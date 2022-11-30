#' Affiliation of assemblages based on phylogenetic turnover 
#'
#' @details The function calculates the degree of affiliation of each community to the region (or evoregion), in which
#'     that community was classified. If used coupled with a analysis of Principal Coordinates of Phylogenetic Structure (PCPS)
#'     to represent the phylogenetic distances the analysis output of the analysis will represent the degree of membership
#'     of assemblages to each evoregion, as described in Maestri and Duarte 2020.
#' 
#' @param phylo.comp.dist A distance matrix indicating the phylogenetic (or taxonomic/functional) distance composition
#'     among assemblages
#' @param groups A character vector indicating the group of each assemblage. This object can be obtained with 
#'     \code{\link{evoregions}} 
#'
#' @return A list with two matrix, one containing affiliation values and the group in which each cell 
#'     is classified and the other containing cell coordinates
#'
#' @references Maestri, R and Duarte L.d.S. (2020). Evoregions: Mapping shifts in phylogenetic turnover across biogeographic regions.
#'     Methods in Ecology and Evolution, 11, 1652-1662.
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' # First run the classification
#' data(akodon.pa.tree) # occurrence data 
#' data(akodon.newick) # phylogenetic tree
#' regions <- evoregions(comm = akodon.pa.tree, phy = akodon.newick)
#' site.region <- regions$Cluster_Evoregions # classification of each community in regions
#' 
#' # Now use the phylogenetic composition (PCPS) and classification obtained from evoregion to calculate the membership
#' axis_sel <- which(regions$PCPS$prop_explainded >= regions$PCPS$tresh_dist) # selecting only significant PCPS axis
#' PCPS_thresh <- regions$PCPS$vectors[, axis_sel] # only significant axis
#' dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean") # distance matrix based on phylogenetic structure
#' calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,groups = regions$Cluster_Evoregions) # affiliation
#' }
#' 
calc_affiliation_evoreg <- function(phylo.comp.dist, groups){
  
  if(class(phylo.comp.dist) != "dist"){
    stop("phylo.comp.dist might be from class dist")
  }
  if(length(groups) != nrow(as.matrix(phylo.comp.dist))){
    stop("Phylogenetic Distance Matrix and group vectors might have the same sites")
  }
  
  n.groups <- length(as.numeric(levels(groups)))
  comm.groups <- groups
  Gs <- lapply(1:n.groups, function(x){
    which(comm.groups == x)
  })
  names(Gs) <- paste("G", 1:n.groups, sep = "")
  dist.matrix <- as.matrix(phylo.comp.dist)
  
  PGall <- lapply(Gs, function(x){
    dist.matrix[x, x]
  })
  
  PGall_similarity <- lapply(PGall, function(x) 1 - x) # distance to similarity
  
  afilliation_by_grp <- lapply(PGall_similarity, function(x){
    afilliation_by_grp <- matrix(NA, nrow(x), 2, dimnames = list(rownames(x), c("afilliation", "group")))
    for (z in 1:nrow(x)) {
      dis <- as.data.frame(x[z,])[-z,]
      afilliation_by_grp[z, 1] <- mean(dis)
    }
    return(afilliation_by_grp)
  })
  list_afilliation_by_grp <- vector(mode = "list", length = length(afilliation_by_grp))
  for(l in 1:length(afilliation_by_grp)){
    afilliation_by_grp_pad <- scales::rescale(afilliation_by_grp[[l]], c(0, 1))
    afilliation_by_grp_pad[, 2] <- l
    list_afilliation_by_grp[[l]] <- afilliation_by_grp_pad
  }
  
  matrix_afilliation <- do.call(rbind, list_afilliation_by_grp)
  matrix_afilliation_org <- matrix_afilliation[match(rownames(dist.matrix), rownames(matrix_afilliation)), ] # organizing the assemblages in the same sequence as PCPS vectors
  return(matrix_afilliation_org)
}