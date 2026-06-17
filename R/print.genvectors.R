#' @rdname genvectors
#' @encoding UTF-8
#' @export
print.genvectors <- function(x, ...){
  cat("$call:\n")
  cat(deparse(x$call), "\n\n")
  if(!is.null(x$haplotypes)){
    cat("$haplotypes:\n")
    print(x$haplotypes, ...)
    cat("\n$individual.per.haplotype:\n")
    print(x$individual.per.haplotype, ...)
  }
  if(!is.null(x$genetic.distances)){
    cat("$genetic.distances:\n")
    print(x$genetic.distances, ...)
    cat("$genetic.per.locality:\n")
    print(x$genetic.per.locality, ...)
  }
  cat("\n$vectors:\n")
  print(x$vectors, ...)
  cat("\n$values:\n")
  print(x$values, ...)
  cat("\n$correlations:\n")
  print(x$correlations, ...)
  cat("\n$P:\n")
  print(x$P, ...)
  cat("\n$scores:\n")
  print(x$scores, ...)
  if(!is.null(x$model)){
    cat("\n$model:\n")
    print(x$model, ...)
    cat("\n$obs.statistic:\n")
    print(x$obs.statistic, ...)
    cat("\n$p.turnover:\n")
    print(x$p.turnover, ...)
    cat("\n$p.divergence:\n")
    print(x$p.divergence, ...)  
  }
  invisible(x)
}