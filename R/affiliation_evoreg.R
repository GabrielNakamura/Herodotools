#' Affiliation values for assemblages according to phylogenetic turnover 
#'
#' @param phylo.comp.dist A distance matrix indicating the phylogenetic (or taxonomic/functional) distance composition
#'     among assemblages
#' @param groups A character vector indicating the group of each assemblage
#'
#' @return A list with two matrix, one containing affiliation values and the group in which each cell 
#'     is classified and the other containing cell coordinates
#'
#' @export
#'
#' @examples
#' 
#' 
affiliation.evoreg <- function(phylo.comp.dist, groups){
  
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