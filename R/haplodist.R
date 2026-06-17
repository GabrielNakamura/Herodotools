#' @title haplodist 
#' 
#' @description Function to extract haplotypes and compute pairwise distances between haplotypes. This function use the functions \code{\link{haplotype}} of package pegas and \code{\link{dist.dna}} of package ape.
#' 
#' @encoding UTF-8
#' @importFrom ape dist.dna
#' @importFrom pegas haplotype
#' @param x A list with the set of DNA sequences (as an object of class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}.
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specify the evolutionary model to be used to compute pairwise distances from DNA sequences (default dist.model = "N").
#' @param ... Additional arguments to the function \code{\link{dist.dna}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes indices that identify each observation sharing the same haplotype.} 
#' \item{individual.per.haplotype}{A matrix with individuals per haplotype.} 
#' \item{haplotype.distances}{A matrix with pairwise distances between haplotypes.} 
#' @seealso \code{\link{HaploVectors}}
#' @examples 
#' data(segv)
#' haplodist(segv$segv.fas)
#' @export
haplodist <- function(x, dist.model = "N", ...){
  res <- list(call = match.call())
  x <- as.list(x)
  n.ind.x <- length(x)
  names.ind.x <- names(x)
  if(is.null(names.ind.x)){
    stop("\n The names are requerid in the object of class DNAbin or haplotype\n")
  }
  hap.x <- pegas::haplotype(x)
  ind.hap.x <- attr(hap.x, "index")
  n.hap <- length(ind.hap.x)
  names.hap.x <- paste("haplotype", attr(hap.x, "dimnames")[[1]], sep = ".")
  attr(hap.x, "dimnames")[[1]] <- names.hap.x
  names(ind.hap.x) <- names.hap.x
  hap.dist <- as.matrix(ape::dist.dna(hap.x, model = dist.model, ...))
  ih <- matrix(0, n.ind.x, n.hap, dimnames = list(names.ind.x, names.hap.x))
  for(i in 1:n.hap){
    ih[ind.hap.x[[i]],i]<-1
  }
  res$haplotypes <- ind.hap.x
  res$individual.per.haplotype <- ih
  res$haplotype.distances <- hap.dist
  return(res)
}