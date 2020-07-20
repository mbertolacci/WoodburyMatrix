#' @name normal-distribution-methods
#' @title Normal distribution methods for \code{SWoodburyMatrix} objects
#' @description
#' Draw samples and compute density functions for the multivariate normal
#' distribution with an \code{SWoodburyMatrix} object as its covariance matrix.
#'
#' @param x A numeric vector or matrix.
#' @param mean Optional mean vector; defaults to zero mean.
#' @param covariance \code{WoodburyMatrix} object.
#' @param log Logical indicating whether to return log of density.
#' @param n Number of samples to return. If \code{n = 1}, returns a vector,
#' otherwise returns an \code{n} by \code{nrow(W)} matrix.
#' @examples
#' library(Matrix)
#' # Trivial example with diagonal covariance matrices
#' W <- WoodburyMatrix(Diagonal(10), Diagonal(10))
#' x <- rwnorm(10, covariance = W)
#' print(dwnorm(x, covariance = W, log = TRUE))
#' @seealso \link{WoodburyMatrix}
NULL

#' @describeIn normal-distribution-methods Compute the density of the
#' distribution
#' @export
dwnorm <- function(x, mean, covariance, log = FALSE) {
  if (missing(mean)) mean <- 0

  n <- if (is.matrix(x)) ncol(x) else length(x)
  output <- -0.5 * (
    n * log(2 * pi)
    + determinant(covariance, logarithm = TRUE)$modulus
    + mahalanobis(x, mean, covariance)
  )
  attributes(output) <- NULL
  if (log) output else exp(output)
}

#' @describeIn normal-distribution-methods Draw samples from the distribution
#' @export
rwnorm <- function(n, mean, covariance) {
  if (missing(mean)) mean <- 0

  output <- sweep(
    as.matrix(
      .sample_Q(n, covariance@A)
      + covariance@X %*% .sample_Q(n, covariance@B)
    ),
    1,
    mean,
    '+'
  )

  if (n == 1) output[, 1] else t(output)
}

.sample_Q <- function(n, Q) {
  z <- matrix(stats::rnorm(n * nrow(Q)), nrow = nrow(Q))

  if (is(Q, 'denseMatrix')) {
    solve(chol(Q), z)
  } else {
    chol_Q <- Cholesky(Q, LDL = FALSE)
    solve(chol_Q, solve(
      chol_Q,
      z,
      system = 'Lt'
    ), system = 'Pt')
  }
}
