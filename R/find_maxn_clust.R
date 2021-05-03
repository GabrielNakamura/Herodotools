#' Estimate the maximum number of groups in DAPC analysis
#'
#' @param x A data.frame or matrix object containing eigenvectors by sites data.
#' @param threshold Scalar. The number of eigenvectors used to perform classification.
#' @param runs Scalar. Number of times classification will be performed.
#' @param method c("kmeans","ward"). See \code{\link{find.clusters}} of adegenet package.
#' @param stat c("BIC", "AIC", or "WSS"). See \code{\link{find.clusters}} of adegenet package.
#' @param criterion c("diffNgroup", "min","goesup", "smoothNgoesup", or "goodfit"). See \code{\link{find.clusters}} of adegenet package.
#' @param max.nclust Value of set of values defining the maximum number of groups to be evaluated.
#' @param subset number of cells used in the analysis. It is particularly important whenever the total number of cells is large (> 1000).
#' @param confidence.level threshold confidence level used to estimate congruence in the classification pattern.
#' @param esp Coordinates for assemblages
#'
#' @return Matrix containing congruence values ranging between 0-1 for each max.nclust value (see Arguments) and confidence level.
#' 
#' @export
#'
#' @examples
find.max.nclust <- function(x, 
                            threshold, 
                            esp,
                            runs = 100, 
                            method = "kmeans", 
                            stat = "BIC", 
                            criterion = "diffNgroup", 
                            max.nclust = c(10, 15, 20, 25, 30), 
                            subset = 100, confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
                            ) 
  {
  group.affinity <- matrix(NA, subset, subset)
  grouping <- matrix(NA, (subset*(subset-1))/2, runs)
  group.sample <- sample(1:nrow(x), size = subset, replace = FALSE)
  congruence <- matrix(NA, length(max.nclust), length(confidence.level), dimnames = list(paste("max.", 
                                                                                               max.nclust, 
                                                                                               "groups", 
                                                                                               sep=" "), 
                                                                                         paste("Confidence.level=", confidence.level, "%",sep="")))
  for (m in 1:length(confidence.level)){
    for (l in 1:length(max.nclust)){
      for (k in 1:runs){
        clust.vec <- adegenet::find.clusters(x[, 1:threshold],
                                           clust = NULL,
                                           choose.n.clust = FALSE,
                                           n.pca = threshold,
                                           method = method,
                                           stat = stat,
                                           n.iter = 1e7,
                                           criterion = criterion,
                                           max.n.clust = max.nclust[[l]]
        )
        rownames(as.data.frame(clust.vec$grp)) == rownames(esp)
        groups <- as.matrix(clust.vec$grp)
        group.subset <- as.matrix(groups[group.sample,])
        tgroups<-as.matrix(t(group.subset))
        for(j in 1:subset){
          for (i in 1:subset){
            group.affinity[i, j] <- group.subset[i, ] == tgroups[ ,j]
          }
        }
        grouping[,k] <- as.vector(as.dist(group.affinity, diag = FALSE, upper=FALSE))
      }
      cor.grouping <- as.matrix(cor(grouping))
      diag(cor.grouping) <- NA
      congruence[l, m] <- sum(ifelse(rowMeans(cor.grouping, na.rm = T) > confidence.level[[m]], 1, 0))/runs
    }
  }  
  congruence
}