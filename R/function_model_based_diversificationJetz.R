#' Tip-based and model-based diversification metrics 
#'
#' @param W Assemblage occurrence matrix, rows are assemblages and columns are species
#' @param tree Phylogenetic tree in newick format
#' @param biogeo Data frame with one column indicating the Ecoregion of each assemblage
#' @param diversification Character indicating the tip-based diversification metric to be calculated. Default is "jetz"
#' @param type Character indicating the type of calculation to be used in the calculation of tip-based metric, "equal-splits" is the default
#' @param ancestral.area One column data frame indicating the Ecoregion of occurrence of each node (rows)
#'
#' @details model-based measure of diversification is calculated as the proportion of total tip-based diversification
#'     metric corresponded only to the branches of a given species that diversified in the current area of occurrence of that species
#'     \deqn{Jetz_total \times (ED_div/ED_total)}
#'
#' @return A list with the following objects:
#' \itemize{
#'     \item{Jetz_per_spp}{Jetz tip-based metric of diversification for each species}
#'      \item{Jetz_species_site}{Matrix containing
#'         the values of Jetz tip-based diversification metric for each assemblage} 
#'     \item{Jetz_harmonic_mean_site}{Harmonic mean of Jetz tip-based
#'         diversification metric}
#'      \item{model_based_Jetz_species_sites}{Matrix with model-based Jetz metric for each species
#'         in each assemblage} 
#'     \item{model_based_Jetz_harmonic_mean_site}{Harmonic mean of model-based Jetz diversification}
#' }
#' 
#' @author Gabriel Nakamura <gabriel.nakamura.souza@@gmail.com> and Arthur V Rodrigues
#'      
#' @importFrom stats cophenetic    
#' 
#' @examples 
#' \dontrun{
#' W_toy <- matrix(c(0, 1, 1, 0, 1, 
#'     1, 0, 1, 0, 1, 
#'     0, 0, 1, 0, 0),nrow= 3,
#'     ncol= 5,dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
#'                           c(paste("s", 1:5, sep=""))))
#' data("toy_treeEx")
#' toy_treeEx                   
#' biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
#' ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
#' plot(toy_treeEx, show.node.label = TRUE)
#' calc_insitu_diversification(W = W_toy,
#'                             tree = toy_treeEx, 
#'                             ancestral.area = ancestral_area_toy,
#'                              biogeo = biogeo_toy)
#' } 
#' 
#' @export
#'
calc_insitu_diversification <- 
  function(W,
           tree,
           ancestral.area, 
           biogeo, 
           diversification = "jetz",
           type = "equal.splits"){
    
    
    if(all(diversification != c("jetz", "freckleton")) == TRUE){
      stop("Tip-based diversification measures must be one of 'jetz' or 'freckleton'")
    }
    
    # Jetz tip-based diversification
    EDtotal <- picante::evol.distinct(tree = tree, type = type)
    Jetz_total <- 1/EDtotal$w
    names(Jetz_total)<- tree$tip.label
    names_spComm <- colnames(W)
    
    # matrix to receive tip-based diversification metrics Jetz
    matrix_XJetz <- matrix(0,
                           nrow = nrow(W),
                           ncol = ncol(W),
                           dimnames = list(rownames(W), colnames(W)))
    
    
    nodes.list <- get_nodes_info_core(W = W, tree = tree, ancestral.area = ancestral.area, biogeo = biogeo) # calculating basic info for db-diversification metrics
    
    if(any(diversification == "jetz")){
      
      ## Jetz local para cada espécie
      l.Jetz.local <- lapply(nodes.list , function(site){
        lapply(1:length(site$nodes_species), function(j){
          # nodes_div <- nodes.list[[3]]$nodes_species[[1]]
          # sp <- names(nodes.list[[3]]$nodes_species)[1]
          nodes_div <- site$nodes_species[[j]]
          sp <- names(site$nodes_species)[j]
          
          if(length(nodes_div) == 1){
            internal.brlen_div <- tree$edge.length[which(tree$edge[,
                                                                   2] %in% nodes_div)] #branch lenghts (in times) for internal branch lengths of the most ancient ancestral
          } else{
            nodes_div <- nodes_div[-1] #internal nodes that form the path from most ancient ancestral that was presented at Ecoregion of local i to species j
            internal.brlen_div <- tree$edge.length[which(tree$edge[,
                                                                   2] %in% nodes_div)] #ages for internal nodes of nodes_div object
          }
          if (length(internal.brlen_div) != 0) { #starts the calculation for ED measure from Redding et al.
            internal.brlen_div <- internal.brlen_div * switch(type, equal.splits = sort(rep(0.5,
                                                                                            length(internal.brlen_div))^c(1:length(internal.brlen_div))),
                                                              fair.proportion = {
                                                                for (j in 1:length(nodes_div)) {
                                                                  sons <- picante::.node.desc(tree, nodes_div[j])
                                                                  n.descendents <- length(sons$tips)
                                                                  if (j == 1) portion <- n.descendents else portion <- c(n.descendents,
                                                                                                                         portion)
                                                                }
                                                                1/portion
                                                              })
            
          }
          
          if(any(is.na(nodes_div))){
            JetzLocal <- 0.00001
          }else{
            ED_div <- sum(internal.brlen_div,
                          tree$edge.length[ape::which.edge(tree, sp)]) #Local ED - modified ED considering only the edges of ancestors inside the biogeo of local i for species j
            
            EDtotal_spp <- EDtotal$w[which(EDtotal$Species == sp)]
            Jetz_total_spp <- Jetz_total[sp] # Diversification calculated according to Jetz for species j
            JetzLocal <- (Jetz_total_spp*(ED_div/EDtotal_spp)) # Modified Jetz to calculate only for local diversification
            
          }
          JetzLocal #Jetz local diversification
          
        })#lapply j == species
      })#lapply site
      
      # populate the matrix
      for(site in 1:length(l.Jetz.local)){
        # site = 1
        for(sp in 1:length(l.Jetz.local[[site]])){
          # sp = 1
          pres <- which(W[site,] >= 1)
          pres <- names_spComm[pres]
          matrix_XJetz[site, pres[sp]] <- l.Jetz.local[[site]][[sp]] # model-base diversification for each species in each site
        }
      }
      
      
      # Summarazing results for each site ---------------------------------------
      
      ## calculating harmonic mean for Jetz
      
      # sum of communities local diversification
      sum_localDiv <- apply(matrix_XJetz,
                            1,
                            function(x) sum(x[which(x != 0 & x != 0.00001)]))
      
      
      #object to receive Jetz diversification values
      matrix_totalDiv_Jetz <- W
      
      # Total jetz values for each community
      for(i in colnames(W)){
        matrix_totalDiv_Jetz[ , i] <- ifelse(matrix_totalDiv_Jetz[ , i] == 1, 
                                             Jetz_total[i],
                                             0) # input Jetz total diversification values in occurence matrix
      }
      
      #sum of communities total diversification
      sum_totalDiv <- apply(matrix_totalDiv_Jetz,
                            1,
                            function(x) sum(x[which(x != 0)]))
      
      numerator_values <- apply(matrix_totalDiv_Jetz, 1, function(x){
        sum(1/x[which(x != 0)])
      })
      
      denom_values <- apply(matrix_totalDiv_Jetz, 1, function(x){
        length(which(x != 0))
      })
      
      #harmonic mean for Jetz total diversification
      JetzTotalComm_harmonic <- denom_values/numerator_values
      
      #harmonic mean for Jetz local diversification
      JetzLocalComm_harmonic <- (JetzTotalComm_harmonic*(sum_localDiv/sum_totalDiv))
      list_res_jetz <- vector(mode = "list", length = 5)
      list_res_jetz[[1]] <- Jetz_total
      list_res_jetz[[2]] <- matrix_totalDiv_Jetz
      list_res_jetz[[3]] <- JetzTotalComm_harmonic # jetz harmonic mean diversification per community
      list_res_jetz[[4]] <- matrix_XJetz # model-based metric per community per species
      list_res_jetz[[5]] <- JetzLocalComm_harmonic # model-based metric harmonic mean 
      names(list_res_jetz) <- c("Jetz_per_spp",
                                "Jetz_species_site",
                                "Jetz_harmonic_mean_site",
                                "model_based_Jetz_species_sites", 
                                "model_based_Jetz_harmonic_mean_site")
      
    }
    
    
    if(any(diversification == "freck")){
      ## Freckleton local
      # modifyed equation 4 from Freckleton et al (2008) considering only the nodes
      # that diversified in ecoregion of local i, plus one is only a correction of the
      # previous step
      
      Freck_total <- numeric(length = length(tree$tip.label))
      names(Freck_total) <- tree$tip.label
      T_freck <- max(cophenetic(tree))/2
      for(l in 1:length(tree$tip.label)){
        nodes_Freck <- picante::.get.nodes(tree, tree$tip.label[l])
        Freck_total[l] <- length(nodes_Freck)/T_freck
      }
      
      matrix_XFreck<- matrix(0,
                             nrow = nrow(W),
                             ncol = ncol(W),
                             dimnames = list(rownames(W), colnames(W)))
      
      for(site in 1:length(nodes.list)){
        for(sp in 1:length(nodes.list[[site]]$nodes_species)){
          pres<- which(W[site,]>=1)
          pres<- names_spComm[pres]
          nodes_div <- nodes.list[[site]]$nodes_species[[sp]]
          matrix_XFreck[site, pres[sp]] <- ifelse(is.na(nodes_div),
                                                  0.00001,
                                                  (length(nodes_div))/T_freck)[1]
        }
      }
      
      FreckLocalComm_mean <- sapply(1:nrow(matrix_XFreck), function(i){
        pres <- matrix_XFreck[i, ][W[i, ] >= 1]
        mean(pres, na.rm = T)
      })
      
      #objet to receive Jetz diversification values
      matrix_totalDiv_Freck <- W
      
      #substitui a matrix de ocorrencia pelos valores de div total do Freck
      for(i in colnames(W)){
        matrix_totalDiv_Freck[, i] <- ifelse(matrix_totalDiv_Freck[, i] >= 1, Freck_total[i], 0)
      }
      
      FreckTotalComm_mean <- sapply(1:nrow(matrix_totalDiv_Freck), function(i){
        pres <- matrix_totalDiv_Freck[i, ][W[i, ] >= 1]
        mean(pres, na.rm = T)
      })
      
      list_res_freckleton <-  data.frame(FreckTotal = FreckTotalComm_mean,
                                         FreckLocal = FreckLocalComm_mean,
                                         FreckLocalProp = FreckLocalComm_mean/
                                           FreckTotalComm_mean,
                                         row.names = row.names(W))
      
    }
    
    
    
    if(all(diversification == c("jetz", "freckleton")) == TRUE){
      list_res$Freckleton <- list_res_freckleton
      list_res$Jetz <- list_res
      return(list_res)
    }
    
    if(diversification == "jetz"){
      return(list_res_jetz)
    }
    
    if(diversification == "freckleton"){
      return(list_res_freckleton)
    }
    
  }