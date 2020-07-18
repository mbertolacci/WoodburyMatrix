#' Instantiate a matrix
#'
#' This is a generic to represent direct instantiation of an implicitly defined
#' matrix. In general this is a bad idea, because it removes the whole purpose
#' of using an implicit representation.
#'
#' @param x Implicit matrix to directly instantiate.
#' @returns The directly instantiated matrix.
#' @seealso \link{WoodburyMatrix}, \linkS4class{WoodburyMatrix}
#' @export
setGeneric('instantiate', function(x) {
  standardGeneric('instantiate')
})

#' @describeIn instantiate Method for general matrices.
#' @export
setMethod(
  'instantiate',
  c(x = 'GWoodburyMatrix'),
  function(x) {
    solve(x@A) + x@U %*% solve(x@B, x@V)
  }
)

#' @describeIn instantiate Method for symmetric matrices.
#' @export
setMethod(
  'instantiate',
  c(x = 'SWoodburyMatrix'),
  function(x) {
    forceSymmetric(solve(x@A) + x@X %*% solve(x@B, t(x@X)))
  }
)
