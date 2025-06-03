#' Computes the arrival ages of each species in the assemblages 
#'
#' @details This function computes the mean arrival age of species in an assemblage based in an ancestral
#'     area reconstruction. For each assemblage we calculate the arrival and establishment of each ancestor of 
#'     present-day species as showed in Van Dijk et al (2021). We consider an arrival event as being the arrival of 
#'     an ancestor and the establishment of this ancestor until the present day. Arrival events that occurred between the 
#'     last speciation event of that lineage and the present time can be assigned a small value or the age corresponding to 
#'     the half of the lenght of the branch linking the ancestor and the present day species
#'
#' @param W Occurrence matrix, rows are assemblages and columns are species
#' @param tree Phylogenetic tree in newick format
#' @param ancestral.area One column data frame, nodes in row and one column containing the occurrence of the nodes
#' @param biogeo One column data frame, assemblages in rows and one column containing the region in which the assemblage is located 
#' @param age.no.ancestor a character string "recent" (default) or "half.edge" indicating how to deal 
#'     with cases in which the most recent aancestor of a tip node are not on the same region. "recent" 
#'     attributes a small age (10e-5), while "half.edge" attributes the age as half the length of the branch
#'     linking the ancestor to the tip.
#' 
#' @return A list of length two. \item{age_arrival}{A matrix with the arrival age of each species in each assemblage} 
#'     \item{mean_age_arrival}{mean age values for each assemblage}
#' 
#' @author Gabriel Nakamura <gabriel.nakamura.souza@@gmail.com> and Arthur Rodrigues
#' 
#' @references Van Dijk, A.; Nakamura G,; Rodrigues, A.V.; Maestri, R. and Duarte, L.d.S. (2021). Imprints of 
#'     tropical niche conservatism and historical dispersal in the radiation of Tyrannidae (Aves:Passeriformes).
#'     Biol. Journ. Linnean Soc., 134, 57-67.
#' 
#' @export
#'
#' @examples
#'  # hypothetical occurrence matrix with species in columns and assemblages in lines
#'  W_toy<- matrix(c(0, 1, 1, 0, 1,
#'                   1, 0, 1, 0, 1, 
#'                   0, 0, 1, 0, 0),nrow= 3,
#'                   ncol= 5,dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
#'                   c(paste("s", 1:5, sep=""))))
#'  
#'  #toy tree
#'  data(toy_treeEx)
#'  
#'  # hypothetical data indicating the ecoregions of each assemblage
#'  biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
#'  
#'  # hypothetical data indicating the ancestral range of each node
#'  ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
#'  
#'  # caculating age of each assemblage
#'  age_assemblages <- calc_age_arrival(W_toy, toy_treeEx, ancestral_area_toy, biogeo_toy)

calc_age_arrival <- function(W, 
                        tree, 
                        ancestral.area, 
                        biogeo,
                        age.no.ancestor = "recent"){
  
  ages <- abs(ape::node.depth.edgelength(phy = tree)
              -max(ape::node.depth.edgelength(tree)))[-c(1:length(tree$tip.label))]
  ages <- data.frame(age = ages)
  s <- length(tree$tip.label)
  
  n.nodes <- tree$Nnode
  rownames(ages) <- paste("N", (s+1):(s+(n.nodes)), sep = "")
  
  age_arrival<-  matrix(0,
                        nrow = nrow(W),
                        ncol = ncol(W),
                        dimnames = list(rownames(W), colnames(W)))
  
  names_spComm <- colnames(W)
  
  nodes.list <- get_nodes_info_core(
    W = W, 
    tree = tree, 
    ancestral.area = ancestral.area, 
    biogeo = biogeo
    )
  
  # site; species; first node (root in the ecoregion)
  
  ## NA means the time of arrival is unknown
  for(site in 1:length(nodes.list)){
    # site = 1
    for(sp in 1:length(nodes.list[[site]]$nodes_species)){
      # sp = 1
      pres<- which(W[site,]>=1)
      pres<- names_spComm[pres]
      node <- nodes.list[[site]]$nodes_species[[sp]][1]
      node.name  <- paste0("N", node)
      age_arrival[site, pres[sp]] <- ages[node.name,] #recebe a idade de chegada
    }
  }
  
  ## age.no.ancestor is an option to atribute an age time when diversification
  ## wasn't local
  if(age.no.ancestor == "recent"){
      age_arrival[is.na(age_arrival)] <- 10e-5
  }
  
  if(age.no.ancestor == "half.edge"){
    for(site in 1:length(nodes.list)){
      for(sp in 1:length(nodes.list[[site]]$disp.anc.node)){
        
        pres<- which(W[site,]==1)
        pres<- names_spComm[pres]
        if(is.na(age_arrival[site,pres[sp]])){
          node <- nodes.list[[site]]$disp.anc.node[sp]
          node.name  <- paste0("N", node)
          age.disp <- ages[node.name,]
          age.arri <- age_arrival[site,pres[sp]] #recebe a idade de chegada
          age.arri <- ifelse(is.na(age.arri), 0, age.arri)
          
          age_arrival[site,pres[sp]] <- mean(c(age.disp, age.arri))
        }
      }
    }
  }
  
  mean_age_arrival <- sapply(1:nrow(age_arrival), function(i){
    pres <- age_arrival[i,][W[i,]==1]
    mean(pres, na.rm = T)
  })
  mean_age_arrival <- data.frame(mean_age_arrival = mean_age_arrival,
                                 row.names = row.names(W))
  list_res <- vector(mode = "list", length = 2)
  list_res[[1]] <- age_arrival
  list_res[[2]] <- mean_age_arrival
  names(list_res) <- c("age_arrival_assemblage", "mean_age_per_assemblage")
  return(list_res)
}
