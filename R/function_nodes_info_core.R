#' Auxiliary function to compute information of node path and dispersal path for each species 
#'
#' @param W Numerical matrix, rows are assemblages and columns are species
#' @param tree Phylogenetic tree in newick format
#' @param biogeo Data frame containing the information of the biome/area/ecoregion of each assemblage
#' @param ancestral.area Single column data frame with nodes in rows and ancestral area/Ecoregion of occurrence in columns. If the reconstruction
#'     comes from BioGeoBears this data frame can be obtained with (\code{\link{get_node_range_BioGeoBEARS}})
#'
#' @return A list with auxiliary information on nodes and species
#' @export
#'
#' @examples
#' \dontrun{
#' # hypothetical occurrence matrix with species in columns and assemblages in lines
#'  W_toy<- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0),
#'      nrow= 3,ncol= 5,
#'      dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
#'      c(paste("s", 1:5, sep=""))))
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
#' get_nodes_info_core(W=W_toy,tree=toy_treeEx,ancestral.area=ancestral_area_toy,biogeo=biogeo_toy)
#' }
#' 
get_nodes_info_core <- function(W,
                            tree,
                            ancestral.area, 
                            biogeo){
  
  names_spComm <- colnames(W)
  AS <- calc_ancestral_state(tree = tree, ancestral.area = ancestral.area)
  nodes.list <- lapply(1:nrow(W), function(i){
    pres <- which(W[i, ] == 1)
    pres <- names_spComm[pres]
    nodes_species <- vector(mode= "list")
    disp.anc.node <- vector("numeric", length = length(pres))
    
    for(j in 1:length(pres)){
      # j = 2
      nodes_sp <- AS[, pres[j]][!is.na(AS[, pres[j]])]
      nodes_sp <- nodes_sp[length(nodes_sp):1] #nodes for species j in community i
      nodes.T <- grep(biogeo[i, 1], nodes_sp)
      nodes_all <- numeric(length = length(nodes_sp))
      nodes_all[nodes.T] <- 1
      
      if(all(nodes_all == 1)){ #if all ancestors of species j are in the same ecoregion of local 1 this will be TRUE
        x <- names(nodes_sp[length(nodes_sp)]) # if TRUE, take basal node as the reference node for calculation of local diversification
        
        rec.anc.node <- as.numeric(substr(x, 2, nchar(x))) # number of ancestral node
        
        n.path <- ape::nodepath(tree,
                                rec.anc.node,
                                which(tree$tip.label == pres[j]))
        
        nodes_species[[j]] <- sort(n.path[-length(n.path)]) # node path from root to tip
        disp.anc.node[j] <- NA # no dispersal
        
      }else
      { # node for the ancestor previous to dispersal
        out.situ.anc.pos <- which(nodes_all == 0)[1]
        x <- names(nodes_sp)[out.situ.anc.pos]
        disp.anc.node[j] <- as.numeric(substr(x, 2, nchar(x)))
        
        # node for the early ancestor in the same ecoregion
        x1 <- names(nodes_sp)[out.situ.anc.pos - 1]
        rec.anc.node <- ifelse(length(x1) == 1,
                               yes = as.numeric(substr(x1, 2, nchar(x1))),
                               no = NA)
        if(!is.na(rec.anc.node)){
          n.path <- ape::nodepath(tree,
                                  rec.anc.node,
                                  which(tree$tip.label==pres[j]))
          
          nodes_species[[j]] <- sort(n.path[-length(n.path)])
        }else{ nodes_species[[j]] <- NA} # options to half length should be implemented here
        
      }
    }
    names(nodes_species) <- pres
    names(disp.anc.node) <- pres
    list(nodes_species = nodes_species,
         disp.anc.node = disp.anc.node)
  })
  return(nodes.list)
}


