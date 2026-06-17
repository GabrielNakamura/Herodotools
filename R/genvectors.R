#' @title GenVectors
#' 
#' @description Function to extract haplotypic/genetic eigenvectors and perform null model-based tests.
#' 
#' @encoding UTF-8
#' @import PCPS
#' @aliases GenVectors print.GenVectors
#' @param pop A matrix describing the incidence of each individual (columns) in a given locality (rows).
#' @param distances Matrix containing genetic distances between individuals or a list with the set of DNA sequences (class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}.
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specify the evolutionary model to be used to computes pairwise distances from DNA sequences (default dist.model = "N").
#' @param log.frequencies Logical argument (TRUE or FALSE) to specify if transformation of natural logarithms plus one in haplotype per locality data must be applied  (Default log.frequencies = FALSE).
#' @param checkdata Logical argument (TRUE or FALSE) to check if individual sequences in the pop data follow the same order as in the set of DNA sequences (Default checkdata = TRUE).
#' @param method Dissimilarity index to apply in matrix P, which describes localities by their haplotypic/genetic composition, as accepted by vegdist function in vegan package (Default method = "euclidean").
#' @param squareroot.dis Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index in matrix P (Default squareroot.dis = TRUE).
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1, 2)).
#' @param analysis Type of analysis, partial match to "none", "adonis" or "glm" (Default analysis = "none"). 
#' @param envir A matrix with environmental variables for each population, with variables as columns and localities as rows. See Details and Examples.
#' @param formula An object of class \code{\link{formula}}. Used in "adonis" or "glm" analysis. See Details and Examples.
#' @param runs Number of permutations for assessing probability of type I error.
#' @param ... Additional arguments to function \code{\link{matrix.p.sig}} and \code{\link{pcps.sig}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes index that identify each observation that share the same haplotype.} 
#' \item{genetic.distances}{A matrix with pairwise genetic/haplotypes distances.}
#' \item{individual.per.haplotype}{A matrix with individuals per haplotype.} 
#' \item{genetic.per.locality}{A matrix with frequency of each genetic/haplotype per locality (\bold{\eqn{W}}).} 
#' \item{vectors}{Haplotypic/genetic eigenvectors.} 
#' \item{values}{Eigenvalues, relative eigenvalues and cumulative relative eigenvalues.} 
#' \item{correlations}{Correlations between haplotypic/genetic eigenvectors and haplotypes/alleles.} 
#' \item{P}{Matrix of haplotypic/genetic composition (\bold{\eqn{P}}).} 
#' \item{scores}{Scores for biplots.} 
#' \item{model}{The observed model.}
#' \item{fun}{The funtion used.}
#' \item{statistic.null.turnover}{A matrix with null statistic for turnover null model.}
#' \item{statistic.null.divergence}{A matrix with null statistic for divergence null model.}
#' \item{statistic.obs}{Observed statistic, F value to predefined function.}
#' \item{p.turnover}{The p value for the turnover null model.}
#' \item{p.divergence}{The p value for the divergence null model.}
#' @details
#' genvectors is a function to extract haplotypic/genetic eigenvectors and perform null model-based 
#' tests. Input parameters can be entered in two ways. The \emph{distances} argument can be supplied 
#' as a set of DNA sequences (class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}, 
#' from which pairwise distances are calculated using \code{\link{dist.dna}}. Alternatively, it can 
#' be provided directly as a genetic distance matrix between individuals.
#' 
#' The argument \emph{analysis} specifies the type of analysis performed. When \emph{analysis} is set
#'  to "adonis" the analysis is performed on a matrix of haplotypic/genetic composition 
#'  (using the \code{\link{matrix.p.sig}} function). The argument \emph{formula} must be specified, 
#'  where the left-hand side gives the resemblance data, right-hand side gives the variables. 
#'  The resemblance data is internally named \emph{p.dist}, thus formula is an expression of the 
#'  form \emph{p.dist ~ predictors}. If \emph{analysis} is set to "glm" it is performed with 
#'  geneticvector (using the \code{\link{pcps.sig}} function). In this case, the argument \emph{formula} 
#'  must also be specified, where the left-hand side gives the vectors used, right-hand side gives the 
#'  variables. The vectors are internally named sequentially \emph{geneticvector.1}, 
#'  \emph{geneticvector.2}, \emph{geneticvector.3}, and so on. Thus, formula is an expression of 
#'  the form \emph{geneticvector.1 ~ predictors}. 
#' 
#' @seealso \code{\link{haplodist}}, \code{\link{matrix.p.sig}}, \code{\link{pcps.sig}}
#' @examples 
#' data(segv)
#' 
#' genvectors(segv$segv.pi, segv$segv.fas, envir = segv$segv.envir,
#'            choices = c(1,2))
#' 
#' genvectors(segv$segv.pi, segv$segv.fas, analysis = "adonis",
#'            envir = segv$segv.envir, formula = p.dist~R, runs = 99)
#' 
#' genvectors(segv$segv.pi, segv$segv.fas, analysis = "glm",
#'            envir = segv$segv.envir, formula = geneticvector.1~R, runs = 99)
#' 
#' 
#' @export
genvectors <- function(pop, distances, checkdata = TRUE, dist.model = "N", log.frequencies = FALSE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", envir, formula, runs = 999, ...){
  if(!inherits(distances, what = "dist") || !inherits(distances, what = "matrix")){
    res <- haplodist(distances, dist.model)
    if (checkdata) {
      if (is.null(colnames(pop))) {
        stop("\n Erro in column names of pop data\n")
      }
      match.names <- match(rownames(res$individual.per.haplotype), colnames(pop))
      if (sum(is.na(match.names)) > 0) {
        print("There are individuals from x data that are not on pop:", quote = FALSE)
        print(setdiff(rownames(res$individual.per.haplotype), colnames(pop)))
        stop("\n Individuals not found on pop matrix\n")
      }
      pop <- as.matrix(pop[, match.names, drop = FALSE])
    } else{
      pop <- as.matrix(pop)
    }
    res$genetic.distances <- res$haplotype.distances
    res$genetic.per.locality <- pop%*%res$individual.per.haplotype
    if(log.frequencies){
      res$genetic.per.locality <- log(res$genetic.per.locality+1)
    }
  } else{
    res <- list()
    if (checkdata) {
      if (is.null(rownames(distances))) {
        stop("\n Erro in row names of distances\n")
      }
      if (is.null(colnames(distances))) {
        stop("\n Erro in column names of distances\n")
      }
      if (is.null(colnames(pop))) {
        stop("\n Erro in row names of pop\n")
      }
      match.names <- match(colnames(pop), colnames(distances))
      if (sum(is.na(match.names)) > 0) {
        stop("\n There are individuals from pop data that are not on distances matrix\n")
      }
      distances <- as.matrix(distances[match.names, match.names])
    }
    res$genetic.distances <- distances
    res$genetic.per.locality <- pop
  }
  res$call <- match.call()
  res.eigen <- PCPS::pcps(res$genetic.per.locality, phylodist = res$genetic.distances, checkdata = FALSE, method = method, squareroot = squareroot.dis)
  res.eigen$call <- NULL
  rownames(res.eigen$values) <- sub("pcps", "geneticvector", rownames(res.eigen$values))
  colnames(res.eigen$vectors) <- sub("pcps", "geneticvector", colnames(res.eigen$vectors))
  colnames(res.eigen$correlations) <- sub("pcps", "geneticvector", colnames(res.eigen$correlations))
  res <- c(res, res.eigen)
  res$scores <- summary(res.eigen, choices = choices)$scores$scores.species
  Analysis <- c("none", "adonis", "glm")
  Analysis <- pmatch(analysis, Analysis)
  if (length(Analysis) != 1 | (is.na(Analysis[1]))) {
    stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
  }
  if(Analysis==2){
    analysis <- "adonis2.margin"
  }
  FUN <- PCPS::select.pcpsmethod(analysis)
  if(Analysis!=1 & !is.null(FUN)){
    if(Analysis == 2){
      test <- PCPS::matrix.p.sig(res$genetic.per.locality, phylodist = res$genetic.distances, checkdata = FALSE, method.p = method, sqrt.p = squareroot.dis, FUN = FUN, envir = envir, runs = runs, newname = "geneticvector", formula = formula, ...)
    }
    if(Analysis == 3){
      choices.analysis <- PCPS::check.formula(formula, colnames(res$vectors))
      test <- PCPS::pcps.sig(res$genetic.per.locality, phylodist = res$genetic.distances, checkdata = FALSE, method = method, squareroot = squareroot.dis, choices = choices.analysis, FUN = FUN, envir = envir, runs = runs, newname = "geneticvector", formula = formula, ...)
    }
    test$call <- NULL
    test$PCPS.obs <- NULL
    names(test) <- sub("null.site", "null.turnover", names(test))
    names(test) <- sub("null.taxa", "null.divergence", names(test))
    names(test) <- sub("site.shuffle", "turnover", names(test))
    names(test) <- sub("taxa.shuffle", "divergence", names(test))
    res <- c(res, test)
  }
  class(res) <- "genvectors"
  return(res)
}