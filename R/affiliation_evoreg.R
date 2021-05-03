#' Affiliation values for assemblages according to phylogenetic turnover 
#'
#' @param evo.vectors An object returned from evoregions function
#' @param n.groups Scalar indicating the number of groups that will be evaluated for affiliation
#' @param method Character indicating the method used to compute the distance in fuzzy matrix
#'
#' @return A list with two elements containing affiliation values and species coordinates
#' @export
#'
#' @examples
afilliation.evoreg <- function(evo.vectors, n.groups = 4, method = "euclidean"){
  groups.vec.bray <- evo.vectors[[2]]$grp
  vec.bray <- evo.vectors[[1]]
  Gs <- lapply(1:n.groups, function(x){
    which(groups.vec.bray == x)
  })
  names(Gs) <- paste("G", 1:n.groups, sep = "")
  
  dist.P.fuzzy <- as.matrix(vegan::vegdist(
                            x = vec.bray,
                            method = "euclidean",
                            diag = T,
                            upper = T
                            )
                            )
    
  PGall <- lapply(Gs, function(x){
    dist.P.fuzzy[x, x]
  })
  fuzzy.OGU.x <- lapply(PGall, function(x){
    fuzzy.OGU.x <- matrix(NA, nrow(x), 2, dimnames = list(rownames(x), c("fuzzy.belonging", "group")))
    for (z in 1:nrow(x)) {
      dis <- as.data.frame(x[z,])[-z,]
      fuzzy.OGU.x[z, 1] <- mean(dis)
    }
    return(fuzzy.OGU.x)
  })
  list_fuzzy.OGU.all.pad <- vector(mode = "list", length = length(fuzzy.OGU.x))
  for(i in 1:length(fuzzy.OGU.x)){
    fuzzy.OGU.x.pad <- scales::rescale(fuzzy.OGU.x[[i]], c(0, 1))
    fuzzy.OGU.x.pad[[i]][, 2] <- i
    list_fuzzy.OGU.all.pad[[i]] <- fuzzy.OGU.x.pad
  }
  fuzzy.OGU.all.pad <- do.call(rbind, list_fuzzy.OGU.all.pad)
  
  org <- SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad, envir = evo.vectors[[3]])
  fuzzy.OGU.all.pad.org <- org$community
  fuzzy.OGU.all.pad.org[, 1] <- 1 - fuzzy.OGU.all.pad.org[, 1]
  esp.fuzzy <- org$environmental
  rownames(fuzzy.OGU.all.pad.org) == rownames(esp.fuzzy)
  list_res <- vector(mode = "list", length = 2)
  list_res[[1]] <- esp.fuzzy
  list_res[[2]] <- fuzzy.OGU.all.pad.org
}