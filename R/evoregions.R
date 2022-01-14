#' Defining groups based on phylogenetic turnover
#'
#' @param comm Species occurrence matrix. Assemblages are rows and species are columns
#' @param coords Numerical matrix with coordinates for assemblages in comm
#' @param phy Newick object with phylogenetic tree containing the species in comm
#' @param method.dist Character. The method to be used to compute phyogenetic distances among assemblages.
#'      Default is "bray"
#' @param tresh.dist A scalar informing the threshold value to select eigenvectors based on the amount of variation of each eigenvector. 
#'     Default is 0.05 (eigenvector with at least 5% of variation)
#' @param method.clust Character indicating the grouping algorithm to be used in cluster analysis. 
#'     Defalt is "kmeans"
#' @param stat.clust Character to be used in \code{\link{find.clusters}} indicating the statistic to be computed for each model. Can be 
#'     "BIC" (default), "AIC" or "WSS". 
#' @param n.iter.clust Integer to be used in \code{\link{find.clusters}} 
#'     function of adegenet package to indicate the number of iterations to be used in each run of K-means algorithm
#' @param criterion.clust Character
#' @param max.n.clust Integer value to be used in \code{\link{find.clusters}}. 
#'     Indicates the maximum number of clusters to be tried. 
#'
#' @return A list of length four containing
#'     1 - Kstat: Numeric vector with the values of the summary statistics for different K values
#'     2 - stat: Numeric value representing the summary statistic for the retained model
#'     3 - grp: A list with three components, the K statistic for each group, the group that each 
#'         assemblage was classified
#'     4 - Size: The number of assemblages in each group
#' 
#' @export
#'
#' @examples
evoregions <- function(comm, 
                       coords, 
                       phy, 
                       method.dist = "bray",
                       tresh.dist = 0.05, 
                       method.clust = "kmeans",
                       stat.clust = "BIC", 
                       n.iter.clust = 1e7, 
                       criterion.clust = "diffNgroup",
                       max.n.clust = 10)
{
  
  if(ape::is.ultrametric(phy) != TRUE){
    stop("Phylogeny must be ultrametric")
  }
  
  zero.comm <- which(rowSums(comm) <= 2) # remove cells with zero, one and two spp
  comm.clean <- comm[-zero.comm, ]
  esp <- coords[-zero.comm, ]
  match <- picante::match.phylo.comm(phy, comm.clean) #standardize species in phylo and comm
  phy <- match$phy
  comm <- match$comm
  zero.comm <- which(rowSums(comm) <= 2) # remove cells with fewer than three species
  comm <- comm[-zero.comm, ]
  esp <- esp[-zero.comm,]
  pcps.comm.bray <- PCPS::pcps(comm, phylodist = cophenetic(phy), method = method.dist)
  P <- pcps.comm.bray$P
  values.bray <- pcps.comm.bray$values
  thresh.bray <- max(which(values.bray[, 2] >= tresh.dist))
  cum.sum.thresh.bray <- cumsum(as.data.frame(values.bray[, 2])
  )[1:thresh.bray, ][3]
  vec.bray <- pcps.comm.bray$vectors
  find.max.number.cluster <- find.max.nclust(x = vec.bray[, 1:thresh.bray],
                                             threshold = thresh.bray, 
                                             nperm = 1000, 
                                             max.nclust = c(10, 15, 20, 25, 30),
                                             subset = 350, 
                                             confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
  )
  clust.vec.bray <- adegenet::find.clusters(vec.bray[, 1:thresh.bray], 
                                            clust = NULL, 
                                            choose.n.clust = FALSE, 
                                            n.pca = thresh.bray, 
                                            method = method.clust, 
                                            stat = stat.clust, 
                                            n.iter = n.iter.clust, 
                                            criterion = criterion.clust, 
                                            max.n.clust = max.n.clust)
  list_res <- vector(mode = "list", length = 3)
  list_res[[1]] <- vec.bray
  list_res[[2]] <- clust.vec.bray
  list_res[[3]] <- esp
  names(list_res) <- c("PCPS_vectors", "Cluster_Evoregions", "Coordinates")
  return(list_res)
}