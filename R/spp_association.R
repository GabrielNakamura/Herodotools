#' Quantifying species association with evoregions
#'
#' @param comm A matrix with assemblages in rows and species in columns
#' @param evo.vectors A list returned from evoregions function
#' @param spp.association A scalar indicating the level of association that will be evaluated. 
#'     Default is 0.7
#'
#' @return A matrix with species association value for each assemblage
#' @export
#'
#' @examples
#' \dontrun{
#' data(akodon.pa.tree) # occurrence data 
#' data(akodon.newick) # phylogenetic tree
#' regions <- evoregions(comm = akodon.pa.tree, phy = akodon.newick) # compute evoregions
#' spp_association_evoreg(comm = akodon.pa.tree, evo.vectors = regions, spp.association = 0.7) # calculate species association
#' }
#' 
spp_association_evoreg <- function(comm, evo.vectors, spp.association = 0.7){
  groups.vec.bray <- evo.vectors[[2]]$grp
  n.groups <- length(evo.vectors[[2]]$size)
  dummy.groups.vec.bray <- matrix(NA, nrow = nrow(groups.vec.bray), ncol = n.groups)
  rownames(dummy.groups.vec.bray) <- rownames(groups.vec.bray)
  colnames(dummy.groups.vec.bray) <- 1:n.groups
  for (i in 1:n.groups) {
    dummy.groups.vec.bray[, i] <- as.numeric(groups.vec.bray == i)
  }
  c <- t(comm) %*% dummy.groups.vec.bray
  c.pad <- vegan::decostand(c, "total")
  spp.groups.bray <- ifelse(c.pad >= spp.association, 1, 0)
  d <- matrix(0, nrow(spp.groups.bray), ncol = 1)
  rownames(d) <- rownames(spp.groups.bray)
  for (k in 1:ncol(spp.groups.bray)) {
    d <- ifelse(spp.groups.bray[, k] == 1, k, d)
  }
  d <- ifelse(d == 0, n.groups + 1, d)
  return(d)
}

