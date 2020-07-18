#' Create a Woodbury matrix identity matrix
#'
#' Creates an implicitly defined matrix representing the equation
#' \deqn{A^{-1} + U B^{-1} V,}
#' where \eqn{A, U, B} and \eqn{V} are n x n, n x p, p x p and p x n matrices,
#' respectively. A symmetric special case is also possible with
#' \deqn{A^{-1} + X B^{-1} X',}
#' where \eqn{X} is n x p and \eqn{A} and \eqn{B} are additionally symmetric.
#' The available methods are described in \link{WoodburyMatrix-class}
#' and in \link{solve}.
#'
#' @details The benefit of using an implicit representation is that the inverse of this
#' matrix can be efficiently calculated via
#' \deqn{A - A U O^{-1} V A}
#' where \eqn{O = B + VAU}, and its determinant by
#' \deqn{det(O) det(A)^{-1} det(B)^{-1}.}
#' These relationships are often called the Woodbury matrix identity and the
#' matrix determinant lemma, respectively. If \eqn{A} and \eqn{B} are sparse or
#' otherwise easy to deal with, and/or when \eqn{p < n}, manipulating the
#' matrices via these relationships rather than forming $W$ directly can have
#' huge advantageous because it avoids having to create the (typically dense)
#' matrix \deqn{A^{-1} + U B^{-1} V} directly.
#'
#' Where applicable, it's worth using the symmetric form of the matrix. This
#' takes advantage of the symmetry where possible to speed up operations, takes
#' less memory, and sometimes has numerical benefits. This function will create
#' the symmetric form in the following circumstances:
#' \itemize{
#'   \item \code{symmetry = TRUE}; or
#'   \item the argument \code{X} is provided; or
#'   \item \code{A} and \code{B} are symmetric (according to
#'   \code{\link{isSymmetric}}) and the arguments \code{U} and \code{V} are
#'   NOT provided.
#' }
#' @returns A \code{\linkS4class{GWoodburyMatrix}} object for a non-symmetric
#' matrix, \code{\linkS4class{SWoodburyMatrix}} for a symmetric matrix.
#' @param A Matrix \eqn{A} in the definition above.
#' @param B Matrix \eqn{B} in the definition above.
#' @param U Matrix \eqn{U} in the definition above. Defaults to a diagonal matrix.
#' @param V Matrix \eqn{V} in the definition above. Defaults to a diagonal matrix.
#' @param X Matrix \eqn{X} in the definition above. Defaults to a diagonal matrix.
#' @param O Optional, precomputed value of \eqn{O}, as defined above. THIS IS
#' NOT CHECKED FOR CORRECTNESS, and this argument is only provided for advanced
#' users who have precomputed the matrix for other purposes.
#' @param symmetric Logical value, whether to create a symmetric or general
#' matrix. See Details section for more information.
#' @seealso \linkS4class{WoodburyMatrix}, \link{solve}, \link{instantiate}
#' @export
WoodburyMatrix <- function(
  A,
  B,
  U = Diagonal(nrow(B)),
  V = Diagonal(ncol(B)),
  X = Diagonal(nrow(B)),
  O,
  symmetric
) {
  if (!missing(X) && (!missing(U) || !missing(V))) {
    stop('Cannot provide X at the same time as U or V')
  }

  if (missing(symmetric)) {
    if (!missing(U) || !missing(V)) {
      symmetric <- FALSE
    } else if (!missing(X)) {
      symmetric <- TRUE
    } else {
      symmetric <- isSymmetric(A) && isSymmetric(B)
    }
  }

  to_Matrix <- function(x) {
    as(x, if (symmetric) 'symmetricMatrix' else 'Matrix')
  }

  A <- to_Matrix(A)
  B <- to_Matrix(B)
  if (!missing(U)) U <- as(U, 'Matrix')
  if (!missing(V)) V <- as(V, 'Matrix')
  if (!missing(X)) X <- as(X, 'Matrix')
  if (missing(O)) {
    if (!symmetric) {
      O <- B + V %*% A %*% U
    } else {
      O <- B + forceSymmetric(crossprod(X, A %*% X))
    }
  }
  O <- to_Matrix(O)

  if (!symmetric) {
    new(
      'GWoodburyMatrix',
      A = A,
      B = B,
      U = U,
      V = V,
      O = O,
      Dim = A@Dim,
      Dimnames = A@Dimnames
    )
  } else {
    new(
      'SWoodburyMatrix',
      A = A,
      B = B,
      X = X,
      O = O,
      Dim = A@Dim,
      Dimnames = A@Dimnames
    )
  }
}
