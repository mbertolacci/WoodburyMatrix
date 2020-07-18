#' @describeIn WoodburyMatrix-class Matrix multiplication (generally fast and
# efficient)
#' @export
setMethod(
  '%*%',
  c(x = 'GWoodburyMatrix', y = 'ANY'),
  function(x, y) {
    solve(x@A, y) + x@U %*% solve(x@B, x@V %*% y)
  }
)

#' @describeIn WoodburyMatrix-class Matrix multiplication (generally fast and
# efficient)
#' @export
setMethod(
  '%*%',
  c(x = 'SWoodburyMatrix', y = 'ANY'),
  function(x, y) {
    solve(x@A, y) + x@X %*% solve(x@B, crossprod(x@X, y))
  }
)

#' @describeIn WoodburyMatrix-class Return the transpose of the matrix as
#' another GWoodburyMatrix.
#' @export
setMethod(
  't',
  c(x = 'GWoodburyMatrix'),
  function(x) {
    WoodburyMatrix(
      A = t(x@A),
      B = t(x@B),
      U = t(x@V),
      V = t(x@U)
    )
  }
)

#' @describeIn WoodburyMatrix-class Does nothing, just returns \code{x}.
#' @export
setMethod(
  't',
  c(x = 'SWoodburyMatrix'),
  function(x) {
    x
  }
)
