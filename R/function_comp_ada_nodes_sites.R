
#' Compute a matrix of node occurrences in each site
#'
#' @param phy a phylogenetic tree
#' @param comm a community matrix 
#'
#' @return a data frame or matrix in dense or long format. Default is dense format (sites x nodes)
#' 
#' @export
#'
#' @examples 
#' 
comp_ada_nodes_sites <-
  function(phy, comm, long = FALSE){
    node_samp_mat <- matrix(NA, 
                            nrow = phy$Nnode,
                            ncol = nrow(comm),
                            dimnames = list(phy$node.label, rownames(comm))
                            )
    names_node <- phy$node.label
    for(i in 1:length(names_node)){
      for(j in 1:nrow(comm)){
        comm_samp <- comm[, which(comm[j, ] == 1)]
        samp_nodes <- picante::prune.sample(samp = comm_samp, phylo = phy)$node.label
        node_samp_mat[i, j] <- ifelse(names_node[i] %in% samp_nodes, 1, 0)
      }
    }
    node_samp_mat <- t(node_samp_mat)
    if(long == TRUE){
      node_samp_mat <- phyloregion::dense2long(node_samp_mat)
      return(node_samp_mat)
    } else{
      return(node_samp_mat)
    }
  }

