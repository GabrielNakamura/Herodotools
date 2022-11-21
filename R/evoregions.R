#' Defining biogeographic regions based on phylogenetic turnover
#' 
#' 
#' @details evoregions performs biogeographical regionalization analysis, 
#'     differently from other methods, evoregion uses a phylogenetic turnover 
#'     metric based on fuzzy sets, therefore accounting for characteristics of 
#'     evolutionary history, e.g tree imbalance, that is not accounted by other
#'     metrics of phylogenetic turnover.
#'
#' @param comm Species occurrence matrix. Assemblages are rows and species are
#'   columns
#' @param phy phylogenetic tree of class `phylo` containing the species in 
#'   `comm`.
#' @param method.dist Character. The method to be used to compute phyogenetic 
#'   distances among assemblages. Dissimilarity index, as accepted by 
#'   \code{\link[vegan]{vegdist}} (Default "bray").
#' @param tresh.dist A scalar informing the threshold value to select 
#'   eigenvectors based on the amount of variation of each eigenvector. 
#'   Default is 0.05 (eigenvector with at least 5% of variation). See details. 
#' @param method.clust Character indicating the grouping algorithm to be used 
#'   in cluster analysis. Options avalible are "kmeans" (default) or "ward". 
#' @param stat.clust Character to be used in 
#'   \code{\link[adegenet]{find.clusters}} indicating the statistic to be 
#'   computed for each model. Can be "BIC" (default), "AIC" or "WSS". 
#' @param n.iter.clust Integer to be used in \code{\link[adegenet]{find.clusters}} 
#'   function of adegenet package to indicate the number of iterations to be
#'   used in each run of K-means algorithm
#' @param criterion.clust a character string matching "diffNgroup" (default),
#'   "min", "goesup", "smoothNgoesup", or "goodfit", indicating the criterion 
#'   for automatic selection of the optimal number of clusters. 
#'   See `criterion` argument in \code{\link[adegenet]{find.clusters} for an explanation
#'   of these procedures.
#' @param max.n.clust Integer value to be used in \code{\link[adegenet]{find.clusters}}. 
#'   Indicates the maximum number of clusters to be tried. 
#'
#' @return A list of length four containing:
#' \itemize{
#'    \item `PCPS_vectors` A matrix with PCPS vectors
#'    \item `Cluster_Evoregion` A vector indicating the region in which each 
#'      assemblage was classified
#' }
#' 
#' @seealso \code{\link{find_max_nclust}} to decide the maximum number of clusters to be 
#'   used
#'   
#' @importFrom stats cophenetic
#' 
#' @example 
#' \dontrun{
#' data(akodon.pa.tree) # occurrence data 
#' data(akodon.newick) # phylogenetic tree
#' regions <- evoregions(comm = akodon.pa.tree, phy = akodon.newick)
#' site.region <- regions$Cluster_Evoregions # classification of each community in regions
#' }
#' 
#' @export
#'
evoregions <- function(comm, 
                       phy, 
                       max.n.clust = NULL,
                       method.dist = "bray",
                       tresh.dist = 0.05, 
                       method.clust = "kmeans",
                       stat.clust = "BIC", 
                       n.iter.clust = 1e7, 
                       criterion.clust = "diffNgroup"
)
{
  
  if(ape::is.ultrametric(phy) != TRUE){
    stop("Phylogeny must be ultrametric")
  }
  
  if(ncol(comm) != length(phy$tip.label)){
    stop("The number of species in the 'comm' and 'phy' do not match. Please use picante::match.phylo.comm() for solve for species matching.")
  }
  
  sp_match_1 <- names(comm) %in% phy$tip.label
  sp_match_2 <- phy$tip.label %in% names(comm)
  
  if(!any(c(sp_match_1, sp_match_2))){
    stop("The names of species in the 'comm' and 'phy' do not match. Please use picante::match.phylo.comm() to solve species matching.")
  }
  
  match <- picante::match.phylo.comm(phy, comm) #standardize species in phylo and comm
  phy <- match$phy
  comm <- match$comm
  
  if(is.null(max.n.clust)){
      matrixP <- SYNCSA::matrix.p(comm = comm, phylodist = cophenetic(phy))
      optimal_matrixP <- phyloregion::optimal_phyloregion(
        x = sqrt(vegan::vegdist(matrixP$matrix.P)), 
        method = "average"
      )
      max.n.clust <- optimal_matrixP$optimal$k
  }
  
  
  pcps.comm.bray <- PCPS::pcps(
    comm,
    phylodist = cophenetic(phy), 
    method = method.dist
  )
  
  P <- pcps.comm.bray$P
  values.bray <- pcps.comm.bray$values
  thresh.bray <- max(which(values.bray[, 2] >= tresh.dist))
  cum.sum.thresh.bray <- cumsum(as.data.frame(values.bray[, 2])
  )[1:thresh.bray, ][3]
  vec.bray <- pcps.comm.bray$vectors
  clust.vec.bray <- adegenet::find.clusters(vec.bray[, 1:thresh.bray], 
                                            clust = NULL, 
                                            choose.n.clust = FALSE, 
                                            n.pca = thresh.bray, 
                                            method = method.clust, 
                                            stat = stat.clust, 
                                            n.iter = n.iter.clust, 
                                            criterion = criterion.clust, 
                                            max.n.clust = max.n.clust)
  list_res <- vector(mode = "list", length = 2)
  
  
  list_res[[1]] <- list(
    vectors = vec.bray,
    prop_explainded = values.bray[,2],
    tresh_dist = tresh.dist
  )
  list_res[[2]] <- clust.vec.bray$grp
  
  
  names(list_res) <- c("PCPS", 
                       "Cluster_Evoregions")
  
  class(list_res) <- "evoregion"
  return(list_res)
}