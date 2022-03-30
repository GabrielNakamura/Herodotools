#' Species and their respective ancestral nodes - internal
#'
#' @param tree phylogenetic tree
#'
#' @return
#'
#' @examples
spp_nodes <- function(tree, node.prefix = "N"){
  tree_base <- phylobase::phylo4(tree)
  n.nodes <- tree$Nnode
  n.spp <- length(tree$tip.label)
  spxnode <- matrix(data = 0, 
                    nrow =  length(tree$tip.label), 
                    ncol =  n.nodes, 
                    dimnames = list(tree$tip.label,
                                    paste(node.prefix, 
                                          (n.spp+1):(n.spp+(n.spp-1)),
                                          sep="")
                                    )
                    )
  for (i in 1:n.spp){
    anc_nodes <- phylobase::ancestors(phy = tree_base, node = i)
    anc_nodes <- anc_nodes - length(tree$tip.label)
    spxnode[i, anc_nodes] <- 1
  }
  t(spxnode)
}