#' Defining groups based on phylogenetic turnover
#'
#' @param comm Community occurrence matrix. Assemblages are rows and species are columns
#' @param coords Numerical matrix with coordinates for assemblages in comm
#' @param phy Newick object with phylogenetic tree containing the species in comm
#' @param method.dist Character with the method to be used to compute phyogenetic distances.
#'      Default is "bray"
#' @param tresh.dist A scalar informing the distance treshold value to select eigenvectors. 
#'     Defalt is 0.05 (eigenvector with at least 5% of variation)
#' @param method.clust Character indicating the grouping algorithm to be used in cluster analysis. 
#'     Defalt is "kmeans"
#' @param stat.clust Character
#' @param n.iter.clust Scalar indicating the number of interections to be used in \code{\link{find.clusters}} function of adegenet package
#' @param criterion.clust Character
#' @param max.n.clust Scalar indicating the maximum number of clusters to be returned
#'
#' @return A adgenet object containing the clusters 
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