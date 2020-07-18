#' Mahalanobis Distance
#'
#' Generic for computing the squared Mahalanobis distance.
#'
#' @param x Numeric vector or matrix
#' @param center Numeric vector representing the mean; if omitted, defaults to
#' zero mean
#' @param cov Covariance matrix
#' @param inverted Whether to treat \code{cov} as a precision matrix; must be
#' \code{FALSE} for \code{SWoodburyMatrix} objects.
#' @param ... Passed to the \code{\link[Matrix]{Cholesky}} function.
#' @seealso \link[base]{mahalanobis}
#' @export
setGeneric('mahalanobis')

#' @describeIn mahalanobis Use the Woodbury matrix identity to compute the
#' squared Mahalanobis distance with the implicit matrix as the covariance.
#' @export
setMethod(
  'mahalanobis',
  signature(x = 'ANY', center = 'ANY', cov = 'SWoodburyMatrix'),
  function(x, center, cov, inverted = FALSE, method = 1, ...) {
    stopifnot(!inverted)

    if (is.matrix(x)) {
      if (!missing(center)) {
        x <- sweep(x, 1, center)
      }
      t_x <- t(x)
      colSums(t_x * solve(cov, t_x))
    } else {
      if (!missing(center)) {
        x <- x - center
      }
      as.numeric(crossprod(x, solve(cov, x)))
    }
  }
)
