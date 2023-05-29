#' Quantifying species association with evoregions
#'
#' @param comm An assemblage/community matrix with communities in rows and species in columns
#' @param group.comm A vector containing the classification of each community in groups
#' @param spp.association A scalar indicating the level of association that will be evaluated to assign species to groups. 
#'     Default is 0.7
#'
#' @return A matrix with species in the rows and one columns indicating the name of the group in which the species 
#'     is associated given the level of species association used in the argument spp.association. If the species fail to reach the threshold 
#'     it will be assigned the "fuzzy" status
#' 
#' @author Leandro Duarte, Renan Maestri and Gabriel Nakamura <gabriel.nakamura.souza@@gmail.com>
#' 
#' @references Maestri, R. and Duarte, L.d.S. (2020). Evoregions: Mapping shifts in phylogenetic turnover across biogeographic regions.
#'     Methods Ecol. Evol., 11, 1652-1662.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data(akodon_sites) # occurrence data 
#' akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
#' data(akodon_newick) # phylogenetic tree
#' spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
#' akodon_pa_tree <- akodon_pa[, spp_in_tree]
#' regions <- calc_evoregions(comm = akodon_pa, phy = akodon_newick) # compute evoregions
#' groups <- regions$cluster_evoregions
#' calc_spp_association_evoreg(comm = akodon_pa_tree, 
#'                        group.comm = groups, 
#'                        spp.association = 0.7) # calculate species association
#' }
#' 
calc_spp_association_evoreg <- 
  function(comm, group.comm, spp.association = 0.7){
    groups.vec.bray <- group.comm
    n.groups <- length(levels(group.comm))
    dummy.groups.vec.bray <- matrix(NA, nrow = length(groups.vec.bray), ncol = as.numeric(n.groups))
    rownames(dummy.groups.vec.bray) <- rownames(comm)
    colnames(dummy.groups.vec.bray) <- paste("group", 1:n.groups, sep = ".")
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
    d <- ifelse(d == 0, "fuzzy", d)
    return(d)
  }

