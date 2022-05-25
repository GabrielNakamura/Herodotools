#' Estimate the maximum number of groups in DAPC analysis
#'
#' @param x A data.frame or matrix object containing eigenvectors by sites data.
#' @param threshold Scalar. The number of eigenvectors used to perform classification.
#' @param nperm Scalar. Number of times classification will be performed.
#' @param method c("kmeans","ward"). See \code{\link{find.clusters}} of adegenet package.
#' @param stat c("BIC", "AIC", or "WSS"). See \code{\link{find.clusters}} of adegenet package.
#' @param criterion c("diffNgroup", "min","goesup", "smoothNgoesup", or "goodfit"). See \code{\link{find.clusters}} of adegenet package.
#' @param max.nclust Value of set of values defining the maximum number of groups to be evaluated.
#' @param subset number of cells used in the analysis. It is particularly important whenever the total number of cells is large (> 1000).
#' @param confidence.level threshold confidence level used to estimate congruence in the classification pattern.
#'
#' @return Matrix containing congruence values ranging between 0-1 for each max.nclust value (see Arguments) and confidence level.
#' 
#' @importFrom stats as.dist cor
#' 
#' @export
#'
#' @examples
#' 

find_max_nclust <- function(x, 
                            threshold, 
                            max.nclust,
                            nperm = 100, 
                            method = "kmeans", 
                            stat = "BIC", 
                            criterion = "diffNgroup",
                            subset = 100,
                            confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
) 
{
  
  if(dim(x) > 1000){
    group.affinity <- matrix(NA, subset, subset)
    grouping <- matrix(NA, (subset*(subset-1))/2, nperm)
    group.sample <- sample(1:nrow(x), size = subset, replace = FALSE)
    congruence <- matrix(NA, length(max.nclust), length(confidence.level), dimnames = list(paste("max.", 
                                                                                                 max.nclust, 
                                                                                                 "groups", 
                                                                                                 sep=" "), 
                                                                                           paste("Confidence.level=", confidence.level, "%",sep="")))
  } else{
    subset <- dim(x)[1]
  }
  vec.perm <- 1:nperm
  group.sample <- sample(1:nrow(x), size = subset, replace = FALSE)
  group_perm <- lapply(max.nclust, function(z){
    lapply(vec.perm, function(y){
      clust_vec <- adegenet::find.clusters(
        x = x[, 1:threshold],
        clust = NULL,
        choose.n.clust = FALSE,
        n.pca = threshold,
        method = "kmeans",
        stat = "BIC",
        n.iter = 1e7,
        criterion = criterion,
        max.n.clust = z)
      groups <- as.matrix(clust_vec$grp)
      group.subset <- as.matrix(groups[group.sample,])
      return(group.subset)
    })
  })
  
  names(group_perm) <- c(paste("group.max", max.nclust, sep =""))
  
  bin_matrix <- lapply(group_perm, function(k){
    do.call(cbind, k)
  })
  
  
  group_affinity <- lapply(bin_matrix, function(x) lapply(apply(x, MARGIN = 2, 
                                                                function(z) lapply(z, function(p) p == z)), function(k) do.call(rbind, k)))
  
  means_corGroup <- lapply(lapply(lapply(group_affinity, function(x){
    lapply(x, function(z){
      as.vector(as.dist(z, diag = FALSE, upper = FALSE))
    })
  }), function(y){
    cor_grouping <- cor(do.call(cbind, y))
    diag(cor_grouping) <- NA
    cor_grouping
  }), function(m) rowMeans(m, na.rm = TRUE))
  
  congruence_matrix <- matrix(unlist(lapply(means_corGroup, function(x){
    lapply(confidence.level, function(y){
      sum(ifelse(x >= y, 1, 0))/nperm
    })
  })), nrow = length(max.nclust), ncol = length(confidence.level), 
  dimnames = list(paste("max.n.clus", 
                        max.nclust, sep = ""), 
                  paste("conf.level", confidence.level, sep = "")))
  
  
  congruence_matrix
}