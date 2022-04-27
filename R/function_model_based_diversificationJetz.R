
function(W,
         tree,
         AS, 
         biogeo, 
         diversification = "jetz",
         type = "equal.splits"){
  
  EDtotal <- picante::evol.distinct(tree = tree, type = type)
  Jetz_total <- 1/EDtotal$w
  names(Jetz_total)<- tree$tip.label
  
  
  # matrix to receive tip-based diversification metrics Jetz
  matrix_XJetz <- matrix(0,
                        nrow = nrow(W),
                        ncol = ncol(W),
                        dimnames = list(rownames(W), colnames(W)))
  
  # matrix to receive the results of local Freckleton metric
  # matrix_XFreck<- matrix(0,
  #                        nrow = nrow(W),
  #                        ncol = ncol(W),
  #                        dimnames = list(rownames(W), colnames(W)))
  
  nodes.list <- nodes_info_core(W = W, tree = tree, AS = AS, biogeo = biogeo) # calculating basic info for db-diversification metrics
  
  if(any(diversification == "jetz")){
    
    ## Jetz local para cada espécie
    l.Jetz.local <- lapply(nodes.list , function(site){
      lapply(1:length(site$nodes_species), function(j){
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
                                                                sons <- .node.desc(tree, nodes_div[j])
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
                        tree$edge.length[ape::which.edge(tree, sp)]) #Local ED - modifyed ED considering only the edges of ancestros inside the biogeo of local i for species j
          
          EDtotal_spp <- EDtotal$w[which(EDtotal$Species==sp)]
          Jetz_total_spp<-Jetz_total[sp] #Diversificatio calculated according to Jetz for species j
          JetzLocal<- (Jetz_total_spp*(ED_div/EDtotal_spp)) #Modified Jetz to calculate only for local diversification
          
        }
        JetzLocal #Jetz local diversification
        
      })#lapply j
    })#lapply site
    
    # populate the matrix
    for(site in 1:length(l.Jetz.local)){
      for(sp in 1:length(l.Jetz.local[[site]])){
        pres<- which(W[site,]>=1)
        pres<- names_spComm[pres]
        matrix_XJetz[site, pres[sp]] <- l.Jetz.local[[site]][[sp]]
      }
    }
    
    
    # Summarazing results for each site ---------------------------------------

    ##calculating harmonic mean for Jetz
    
    #sum of communities local diversification
    sum_localDiv <- apply(matrix_XJetz,
                         1,
                         function(x) sum(x[which(x != 0 & x != 0.00001)]))
    
    
    #objet to receive Jetz diversification values
    matrix_totalDiv_Jetz <- W
    
    #substitui a matrix de ocorrencia pelos valores de div total do Jetz
    for(i in colnames(W)){
      #i=1
      matrix_totalDiv_Jetz[ , i] <- ifelse(matrix_totalDiv_Jetz[ , i] >= 1, 
                                           Jetz_total[i],
                                           0) #input Jetz total diversification values in occurence matrix
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
    
  }
  
 #if{any(diversification) == "Freckleton"}{
 #  
 #}
  
  list_res <- vector(mode = "list", length = 4)
  list_res[[1]] <- Jetz_total
  list_res[[2]] <- matrix_XJetz
  list_res[[3]] <- JetzLocalComm_harmonic # model-based metric
  names(list_res) <- c("Jetz_per_spp", "model-basedJetz_species_site", "model-based_totalJetz_persite")
  return(list_res)
}
