#' Compute node ages
#'
#' @param phy a phylogenetic tree with nodes named
#'
#' @return a data frame with two columns, one with species and node names, another with node ages
#'     species will show age 0
#'     
#' @export
#'
#' @examples
comp_age <- function(phy){
  root.age <- max(cophenetic(phy))/2
  ages <- data.frame(names = c(phy$tip.label, phy$node.label),
                     ages = round(root.age - ape::node.depth.edgelength(phy), 5))
  return(ages)
}