#' Defining groups based on phylogenetic turnover
#'
#' @param comm Species occurrence matrix. Assemblages are rows and species are columns
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
#' @param criterion.clust Character indicating the type of method to be used to define the maximum number of clusters. Default option
#'     is "elbow" method, implemented in \code{\link{find.clusters}} of phyloregion package.
#' @param max.n.clust Integer value to be used in \code{\link{find.clusters}}. 
#'     Indicates the maximum number of clusters to be tried. 
#' @param max.n.clust.method 
#'
#' @return A list of length four containing
#'     \item{PCPS_vectors}: A matrix with PCPS vectors
#'     \item{Cluster_Evoregion}: A vector indicating the group in which each assemblage was classified
#'     \item{Coordinates}: A matrix with spatial coordinates of assemblages
#'     \item{Matrix_P}: Phylogenetic composition matrix obtained from [matrix.p()] function from SYNCSA package
#'     \item{Total_vec_var}: A numeric vector with the amount of variation in the PCPS vectors selected accordingly with the threshold argument
#' 
#' @seealso [find.max.nclust()] to decide the maximum number of clusters to be used
#' 
#' @export
#'
#' @examples
#' 

evoregions <- function(comm, 
                       phy, 
                       max.n.clust = NULL,
                       max.n.clust.method = "elbow",
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
  
  zero.comm <- which(rowSums(comm) <= 1) # remove cells with zero, one and two spp
  comm.clean <- comm[-zero.comm, ]
  match <- picante::match.phylo.comm(phy, comm.clean) #standardize species in phylo and comm
  phy <- match$phy
  comm <- match$comm
  if(is.null(max.n.clust)){
    match_clust <- pmatch(max.n.clust.method, "elbow")
    if(is.na(match_clus)){
      stop("wrong method to define the number of clusters")
    }
    if(match_clust == 1){ # elbow method
      matrixP <- SYNCSA::matrix.p(comm = comm, phylodist = cophenetic(phy))
      optimal_matrixP <- phyloregion::optimal_phyloregion(x = sqrt(vegan::vegdist(matrixP$matrix.P)), method = "average")
      max.n.clust <- optimal_matrixP$optimal$k
    }
  }
  pcps.comm.bray <- PCPS::pcps(comm, phylodist = cophenetic(phy), method = method.dist)
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
  list_res <- vector(mode = "list", length = 4)
  list_res[[1]] <- vec.bray
  list_res[[2]] <- clust.vec.bray
  list_res[[3]] <- P
  list_res[[4]] <- cum.sum.thresh.bray
  names(list_res) <- c("PCPS_vectors", "Cluster_Evoregions", "Matrix_P", "Total_vec_var")
  return(list_res)
}