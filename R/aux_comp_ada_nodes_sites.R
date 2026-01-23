#' Axiliary function to compute node composition by sites
#'
#' @param phy A phylo object
#' @param comm A community matrix, rows are communities and columns are species
#' @param long Logical, if TRUE dense format of the output matrix is converted
#'     to long format
#'
#' @returns A matrix (or data frame in the long format) with community composition
#'     by nodes
#' @export

comp_ada_nodes_sites <-
  function(phy, comm, long = FALSE) {
    node_samp_mat <- matrix(
      NA,
      nrow = phy$Nnode,
      ncol = nrow(comm),
      dimnames = list(phy$node.label, rownames(comm))
    )
    
    names_node <- phy$node.label
    
    for (i in 1:length(names_node)) {
      for (j in 1:nrow(comm)) {
        comm_samp <- comm[, which(comm[j, ] == 1), drop = FALSE]
        samp_nodes <- picante::prune.sample(samp = comm_samp, phylo = phy)$node.label
        node_samp_mat[i, j] <- ifelse(names_node[i] %in% samp_nodes, 1, 0)
      }
    }
    
    node_samp_mat <- t(node_samp_mat)
    
    if (long == TRUE) {
      node_samp_mat <- phyloregion::dense2long(node_samp_mat)
      return(node_samp_mat)
    } else {
      return(node_samp_mat)
    }
  }

