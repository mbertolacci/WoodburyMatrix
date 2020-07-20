#' @describeIn WoodburyMatrix-class Check for symmetry of matrix; always returns
#' \code{FALSE}.
#' @export
setMethod(
  'isSymmetric',
  c(object = 'GWoodburyMatrix'),
  function(object) FALSE
)

#' @describeIn WoodburyMatrix-class Check for symmetry of matrix; always returns
#' \code{TRUE}.
#' @export
setMethod(
  'isSymmetric',
  c(object = 'SWoodburyMatrix'),
  function(object) TRUE
)
