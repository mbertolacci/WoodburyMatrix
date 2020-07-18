#' Mahalanobis Distance
#'
#' Generic for computing the squared Mahalanobis distance.
#'
#' @param x Numeric vector
#' @param center Numeric vector representing the mean; if omitted, defaults to
#' zero mean
#' @param cov Covariance matrix
#' @param inverted Whether to treat \code{cov} as a precision matrix; must be
#' \code{FALSE} for \code{SWoodburyMatrix} objects.
#' @param ... Passed to the \code{\link[Matrix]{Cholesky}} function.
#' @seealso base::mahalanobis
#' @export
setGeneric('mahalanobis')

#' @describeIn mahalanobis Use the Woodbury matrix identity to compute the
#' squared Mahalanobis distance with the implicit matrix as the covariance.
#' @export
setMethod(
  'mahalanobis',
  signature(x = 'ANY', center = 'ANY', cov = 'SWoodburyMatrix'),
  function(x, center, cov, inverted = FALSE, ...) {
    stopifnot(!inverted)

    if (!missing(center)) {
      x <- x - center
    }

    chol_O <- .try_Cholesky(cov@O, LDL = FALSE, ...)
    if (!is.null(chol_O)) {
      A_x <- cov@A %*% x
      as.numeric(crossprod(
        x,
        A_x
      ) - crossprod(.chol_Linv_P(
        chol_O,
        crossprod(cov@X, A_x)
      )))
    } else {
      as.numeric(crossprod(x, solve(cov, x)))
    }
  }
)
