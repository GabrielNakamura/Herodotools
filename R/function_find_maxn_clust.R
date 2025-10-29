# auxiliar funciton for computing affinity for each max clust level

calc_pairwise_group_classification <- 
  function(group_matrix){
    lapply(seq_len(ncol(group_matrix)), function(j) {
      as.dist(outer(group_matrix[, j], group_matrix[, j], `==`))
    })
  }


#' Estimate the maximum number of groups in DAPC analysis
#'
#' @details
#' Additional details...
#' 
#' 
#' @param x A data.frame or matrix object containing eigenvectors by sites.
#' @param threshold Scalar. The number of eigenvectors used to perform classification.
#' @param nperm Scalar. Number of times classification will be performed.
#' @param method Character, one of c("kmeans","ward"). This will be used in `find.clusters` function. 
#'     See \code{\link{find.clusters}} of adegenet package. Default is "kmeans"
#' @param stat Character, one of c("BIC", "AIC", or "WSS"). This will be used in `find.clusters` function. 
#'     See \code{\link{find.clusters}} of adegenet package. Default is "BIC".
#' @param criterion Character one of c("diffNgroup", "min","goesup", "smoothNgoesup", or "goodfit"). 
#'     This will be used in `find.clusters` function. Default is "diffNgroup". 
#'     See \code{\link{find.clusters}} of adegenet package.
#' @param max.nclust A vector containing values of the maximum number of groups to be evaluated.
#' @param subset Scalar. The number of cells used in the analysis. 
#'     It is particularly important whenever the total number of cells is large (> 1000).
#' @param confidence.level A vector containing values with threshold confidence
#'     level used to estimate congruence in the classification pattern.
#'
#' @return Matrix containing congruence values ranging between 0-1 for each max.nclust 
#'     value (see Arguments) and confidence level.
#' 
#' @importFrom stats as.dist cor
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data(regions)
#' evovectors <- regions$PCPS$vectors # eigenvectors by site
#' find_max_nclust(x = evovectors, threshold = 3, max.nclust = 10)
#' }

find_max_nclust <- 
  function(x, 
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
    # in case of big matrices, use only a subset of that matrix
    if(dim(x)[1] > 1000){
      group.affinity <- matrix(NA, subset, subset)
      grouping <- matrix(NA, (subset*(subset-1))/2, nperm)
      group.sample <- sample(1:nrow(x), size = subset, replace = FALSE)
      congruence <- matrix(NA, length(max.nclust), length(confidence.level), dimnames = list(paste("max.", 
                                                                                                   max.nclust, 
                                                                                                   "groups", 
                                                                                                   sep=" "), 
                                                                                             paste("Confidence.level=", confidence.level, "%",sep="")))
    } else{ # the whole matrix is used if it is lower than 1000 assemblages
      subset <- dim(x)[1]
    }
    vec.perm <- 1:nperm
    group.sample <- sample(1:nrow(x), size = subset, replace = FALSE)
    
    total_steps <- length(max.nclust) * nperm  # total iterations of inner loop
    
    progressr::with_progress({
      p <- progressr::progressor(steps = total_steps)
      group_perm <- 
        future.apply::future_lapply(max.nclust, function(z) {
          future.apply::future_lapply(vec.perm, function(y) {
            
            # update progress
            p(message = sprintf("max.nclust=%s perm=%s", z, y))
            
            clust_vec <- adegenet::find.clusters(
              x = x[, 1:threshold],
              clust = NULL,
              choose.n.clust = FALSE,
              n.pca = threshold,
              method = method,
              stat = stat,
              n.iter = 1e5,        # Reduced from 1e7 → MUCH faster, same result quality
              criterion = criterion,
              max.n.clust = z      # Now actually uses z
            )
            
            # just in case the group sample is used, if all the samples they are only in a different order
            groups <- clust_vec$grp[group.sample]
            return(groups)
          })
        })
    })
    
    # renaming all classifications for all levels of max cluster and all repetitions
    names(group_perm) <- c(paste("group.max", max.nclust, sep =""))
    
    # binding in a data frame all classifications, each for max level cluster group
    bin_matrix <- 
      lapply(group_perm, function(k){
        do.call(cbind, k)
      })
    
    
    # list with max clust length and for each nruns vectors indicating the affinity (pairwise grouping)
    list_affinity_runs <- 
      lapply(bin_matrix, 
             function(x) as.vector(calc_pairwise_group_classification(x)))
    
    # joining affinities for each level of max clust in a matrix 
    list_matrix_affinities <- 
      lapply(list_affinity_runs, 
             function(x) do.call(cbind, x))
    
    # correlation for each level of max clust
    list_cor_affinities <- 
      lapply(list_matrix_affinities, function(x){
        mat1 <- as.matrix(cor(x))
        diag(mat1) <- NA
        return(mat1)
      })
    
    mean_correlation_group <- lapply(list_cor_affinities, function(x) rowMeans(x, na.rm = T))
    
    matrix_mean_correlation <- do.call(cbind, mean_correlation_group)
    
    list_confidence_test <- 
      apply(matrix_mean_correlation, 2, function(x){
        lapply(confidence.level, function(z) sum(ifelse(x >= z, 1, 0))/nperm )
      })
    
    mat_cor_group_level <- do.call(rbind, lapply(list_confidence_test, function(x) do.call(cbind, x)))
    rownames(mat_cor_group_level) <- names(list_confidence_test)
    colnames(mat_cor_group_level) <- paste("confidence_lev", confidence.level, sep = "_")
    congruence <- mat_cor_group_level
    return(congruence)
    
  }