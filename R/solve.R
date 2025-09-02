#' @name solve-methods
#' @aliases solve
#' @title Solve methods for \code{WoodburyMatrix} objects
#' @description
#' Methods based on \code{\link{solve}} to solve a linear system of equations
#' involving \code{\linkS4class{WoodburyMatrix}} objects. These methods take
#' advantage of the Woodbury matrix identity and therefore can be much more
#' time and memory efficient than forming the matrix directly.
#'
#' Calling this function while omitting the \code{b} argument returns the
#' inverse of \code{a}. This is NOT recommended, since it removes any benefit
#' from using an implicit representation of \code{a}.
#'
#' @param a \code{WoodburyMatrix} object.
#' @param b Matrix, vector, or similar (needs to be compatible with the
#' submatrices \code{a@A} and \code{a@V} or \code{a@X} that define the
#' \code{WoodburyMatrix}).
#' @returns The solution to the linear system, or the inverse of the matrix. The
#' class of the return value will be a vector if \code{b} is a vector, and may
#' otherwise be either a regular matrix or a subclass of
#' \code{\link[Matrix:Matrix-class]{Matrix}}, with the specific subclass determined
#' by \code{a} and \code{b}.
#' @seealso \link{WoodburyMatrix}, \linkS4class{WoodburyMatrix}
NULL

setClassUnion('numLike', members = c('logical', 'numeric'))

#' @describeIn solve-methods Invert the matrix
#' @export
setMethod(
  'solve',
  signature(a = 'GWoodburyMatrix', b = 'missing'),
  function(a) {
    a@A - a@A %*% a@U %*% solve(a@O, a@V %*% a@A)
  }
)

.solve_b_base <- function(a, b) {
  A_b <- a@A %*% b
  A_b - a@A %*% a@U %*% solve(a@O, a@V %*% A_b)
}

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'GWoodburyMatrix', b = 'numLike'),
  function(a, b) {
    drop(.solve_b_base(a, b))
  }
)

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'GWoodburyMatrix', b = 'matrix'),
  function(a, b) {
    .solve_b_base(a, b)
  }
)

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'GWoodburyMatrix', b = 'ANY'),
  function(a, b) {
    .solve_b_base(a, b)
  }
)

#' @describeIn solve-methods Invert the symmetric matrix
#' @export
setMethod(
  'solve',
  signature(a = 'SWoodburyMatrix', b = 'missing'),
  function(a) {
    Xt_A <- crossprod(a@X, a@A)
    a@A - forceSymmetric(crossprod(Xt_A, solve(a@O, Xt_A)))
  }
)

.solve_b_base_symmetric <- function(a, b) {
  A_b <- a@A %*% b
  A_b - a@A %*% a@X %*% solve(a@O, crossprod(a@X, A_b))
}

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'SWoodburyMatrix', b = 'numLike'),
  function(a, b) {
    drop(.solve_b_base_symmetric(a, b))
  }
)

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'SWoodburyMatrix', b = 'matrix'),
  function(a, b) {
    .solve_b_base_symmetric(a, b)
  }
)

#' @describeIn solve-methods Solve the linear system
#' @export
setMethod(
  'solve',
  signature(a = 'SWoodburyMatrix', b = 'ANY'),
  function(a, b) {
    .solve_b_base_symmetric(a, b)
  }
)
