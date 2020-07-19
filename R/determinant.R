#' Calculate the determinant of a WoodburyMatrix object
#'
#' Calculates the (log) determinant of a \code{\linkS4class{WoodburyMatrix}}
#' using the matrix determinant lemma.
#'
#' @param x A object that is a subclass of \code{WoodburyMatrix}
#' @param logarithm Logical indicating whether to return the logarithm of the
#' matrix.
#' @returns Same as \link[base]{determinant}.
#' @seealso \link[base]{determinant}
#' @export
setMethod(
  'determinant',
  signature(x = 'WoodburyMatrix', logarithm = 'logical'),
  function(x, logarithm) {
    A_det <- determinant(x@A, logarithm = TRUE)
    B_det <- determinant(x@B, logarithm = TRUE)
    O_det <- determinant(x@O, logarithm = TRUE)

    modulus <- O_det$modulus - A_det$modulus - B_det$modulus
    sign <- A_det$sign * B_det$sign * O_det$sign

    if (!logarithm) modulus <- exp(modulus)
    attr(modulus, 'logarithm') <- logarithm
    structure(list(modulus = modulus, sign = sign), class = 'det')
  }
)
