#' Calculate the Harmonic Mean
#'
#' Computes the harmonic mean of a numeric vector.  
#' The harmonic mean is defined as \eqn{n / \sum (1/x_i)},  
#' where \eqn{n} is the number of observations.  
#'
#' @param x A numeric vector.
#' @param na.rm Logical; if \code{TRUE}, remove \code{NA} values before computation.
#' @param ignore.zero Logical; if \code{TRUE}, remove zeros before computation.
#'
#' @details
#' - If \code{ignore.zero = FALSE} and the vector contains zeros,
#'   the harmonic mean is zero by definition.
#' - If all values are removed (due to \code{NA}s or zeros),
#'   the function returns \code{NA}.
#'
#' @return A numeric value representing the harmonic mean.
#'
#' @examples
#' # Simple harmonic mean
#' calc_harmonic_mean(c(1, 2, 3))
#'
#' # With a zero (harmonic mean collapses to 0)
#' calc_harmonic_mean(c(1, 2, 0))
#'
#' # Ignoring zeros
#' calc_harmonic_mean(c(1, 2, 0), ignore.zero = TRUE)
#'
#' # With NA values
#' calc_harmonic_mean(c(1, 2, NA), na.rm = TRUE)
#'
#' # Combining both options
#' calc_harmonic_mean(c(1, 2, 0, NA), na.rm = TRUE, ignore.zero = TRUE)
#'
#' @export
calc_harmonic_mean <- function(x, na.rm = FALSE, ignore.zero = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  if (ignore.zero) {
    x <- x[x != 0]
  }
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  length(x) / sum(1 / x)
}
