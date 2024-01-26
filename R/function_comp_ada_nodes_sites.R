get_spp_nodes2 <- 
  function(tree){
    tree_base <- phylobase::phylo4(tree)
    n.nodes <- tree$Nnode
    n.spp <- length(tree$tip.label)
    spxnode <- matrix(data = 0, 
                      nrow =  length(tree$tip.label), 
                      ncol =  n.nodes, 
                      dimnames = list(tree$tip.label,
                                      tree$node.label)
    )
    for (i in 1:n.spp){
      anc_nodes <- phylobase::ancestors(phy = tree_base, node = i)
      anc_nodes <- anc_nodes - length(tree$tip.label)
      spxnode[i, anc_nodes] <- 1
    }
    t(spxnode)
  }

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
    node_comp <- suppressWarnings(get_spp_nodes2(phy))  
    for(j in 1:nrow(comm)){
      comm_samp <- comm[j, which(comm[j, ] == 1)]
      tree_comm <- ape::keep.tip(phy = phy, tip = names(comm_samp))
      nodes_comm <- rownames(get_spp_nodes2(tree_comm))
      node_samp_mat[nodes_comm, j] <- 1
    }
    node_samp_mat <- t(node_samp_mat)
    node_samp_mat <- ifelse(is.na(node_samp_mat), 0, node_samp_mat)
    if(long == TRUE){
      node_samp_mat <- phyloregion::dense2long(node_samp_mat)
      return(node_samp_mat)
    } else{
      return(node_samp_mat)
    }
  }
