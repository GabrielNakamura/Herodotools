evo.vectors <- phylo_evoregion
method <- "euclidean"
dist.matrix <- pb$phylo.beta.sim
groups <- phylo_regionalization$region.df
afilliation.evoreg <- function(evo.vectors, dist.matrix = NULL, phyloregion.classification = NULL,
                               method = "euclidean"){
  if(!is.null(dist.matrix) == TRUE){
    if(!is.null(phyloregion.classification) == TRUE){
      groups.vec.bray <- phylo_regionalization$region.df$cluster
      names(groups.vec.bray) <- phylo_regionalization$region.df$grids
      groups.vec.bray <- as.factor(groups.vec.bray)
      n.groups <- length(levels(as.factor(groups.vec.bray)))
    }
    
  }
  n.groups <- length(as.numeric(levels(evo.vectors[[2]]$grp)))
  groups.vec.bray <- evo.vectors[[2]]$grp
  vec.bray <- evo.vectors[[1]]
  Gs <- lapply(1:n.groups, function(x){
    which(groups.vec.bray == x)
  })
  names(Gs) <- paste("G", 1:n.groups, sep = "")
  dist.P.fuzzy <- as.matrix(vegan::vegdist(
    x = vec.bray,
    method = method,
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
    fuzzy.OGU.x.pad[, 2] <- i
    list_fuzzy.OGU.all.pad[[i]] <- fuzzy.OGU.x.pad
  }
  fuzzy.OGU.all.pad <- do.call(rbind, list_fuzzy.OGU.all.pad)
  
  org <- SYNCSA::organize.syncsa(comm = fuzzy.OGU.all.pad, envir = evo.vectors[[3]])
  fuzzy.OGU.all.pad.org <- org$community
  fuzzy.OGU.all.pad.org[, 1] <- 1 - fuzzy.OGU.all.pad.org[, 1]
  esp.fuzzy <- org$environmental
  list_res <- vector(mode = "list", length = 2)
  list_res[[1]] <- esp.fuzzy
  list_res[[2]] <- fuzzy.OGU.all.pad.org
  list_res
}