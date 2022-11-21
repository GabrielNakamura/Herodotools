#' Species and their respective ancestral nodes
#'
#' @details This is an auxiliary function that, by taking a phylogenetic tree, allows to represent in a matrix format 
#'     the relationship between present-day species (tips) and their ancestors (nodes). The function provides an
#'     easy way to obtain certain objects that can be used in other functions, like \code{\link{ancestral_state}}
#'
#' @description This function computes a matrix containing species and their respective ancestral nodes
#'
#' @param tree Phylogenetic tree
#' @param node.prefix Character indicating the prefix to be used to name nodes, default is the letter "N"
#'
#' @return A matrix with species in lines and nodes in columns. 1 indicates that the node corresponds to the ancestor of that species and zero indicates that species 
#'     does not share the ancestor node
#'     
#' @export
#' 
#' @examples
#' data(akodon.newick)
#' spp_nodes(tree = akodon.newick, node.prefix = "N")
#' 
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