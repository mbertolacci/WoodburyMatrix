#' Create a Woodbury matrix identity matrix
#'
#' Creates an implicitly defined matrix representing the equation
#' \deqn{A^{-1} + U B^{-1} V,}
#' where \eqn{A, U, B} and \eqn{V} are n x n, n x p, p x p and p x n matrices,
#' respectively. A symmetric special case is also possible with
#' \deqn{A^{-1} + X B^{-1} X',}
#' where \eqn{X} is n x p and \eqn{A} and \eqn{B} are additionally symmetric.
#' The available methods are described in \link{WoodburyMatrix-class}
#' and in \link{solve}. Multiple B / U / V / X matrices are also supported; see
#' below
#'
#' @details The benefit of using an implicit representation is that the inverse
#' of this matrix can be efficiently calculated via
#' \deqn{A - A U O^{-1} V A}
#' where \eqn{O = B + VAU}, and its determinant by
#' \deqn{det(O) det(A)^{-1} det(B)^{-1}.}
#' These relationships are often called the Woodbury matrix identity and the
#' matrix determinant lemma, respectively. If \eqn{A} and \eqn{B} are sparse or
#' otherwise easy to deal with, and/or when \eqn{p < n}, manipulating the
#' matrices via these relationships rather than forming \eqn{W} directly can
#' have huge advantageous because it avoids having to create the (typically
#' dense) matrix \deqn{A^{-1} + U B^{-1} V} directly.
#'
#' @section Symmetric form:
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
#'
#' @section Multiple B matrices:
#' A more general form allows for multiple B matrices:
#' \deqn{A^{-1} + \sum_{i = 1}^n U_i B_i^{-1} V_i,}
#' and analogously for the symmetric form. You can use this form by providing
#' a list of matrices as the \code{B} (or \code{U}, \code{V} or \code{X})
#' arguments. Internally, this is implemented by converting to the standard form
#' by letting \code{B = bdiag(...the B matrices...)},
#' \code{U = cbind(..the U matrices...)}, and so on.
#'
#' The \code{B}, \code{U}, \code{V} and \code{X} values are recycled to the
#' length of the longest list, so you can, for instance, provide multiple B
#' matrices but only one U matrix (and vice-versa).
#'
#' @returns A \code{\linkS4class{GWoodburyMatrix}} object for a non-symmetric
#' matrix, \code{\linkS4class{SWoodburyMatrix}} for a symmetric matrix.
#' @param A Matrix \eqn{A} in the definition above.
#' @param B Matrix \eqn{B} in the definition above, or list of matrices.
#' @param U Matrix \eqn{U} in the definition above, or list of matrices.
#' Defaults to a diagonal matrix/matrices.
#' @param V Matrix \eqn{V} in the definition above, or list of matrices.
#' Defaults to a diagonal matrix/matrices.
#' @param X Matrix \eqn{X} in the definition above, or list of matrices.
#' Defaults to a diagonal matrix/matrices.
#' @param O Optional, precomputed value of \eqn{O}, as defined above. THIS IS
#' NOT CHECKED FOR CORRECTNESS, and this argument is only provided for advanced
#' users who have precomputed the matrix for other purposes.
#' @param symmetric Logical value, whether to create a symmetric or general
#' matrix. See Details section for more information.
#' @examples
#' library(Matrix)
#' # Example solving a linear system with general matrices
#' A <- Diagonal(100)
#' B <- rsparsematrix(100, 100, 0.5)
#' W <- WoodburyMatrix(A, B)
#' str(solve(W, rnorm(100)))
#'
#' # Calculating the determinant of a symmetric system
#' A <- Diagonal(100)
#' B <- rsparsematrix(100, 100, 0.5, symmetric = TRUE)
#' W <- WoodburyMatrix(A, B, symmetric = TRUE)
#' print(determinant(W))
#'
#' # Having a lower rank B matrix and an X matrix
#' A <- Diagonal(100)
#' B <- rsparsematrix(10, 10, 1, symmetric = TRUE)
#' X <- rsparsematrix(100, 10, 1)
#' W <- WoodburyMatrix(A, B, X = X)
#' str(solve(W, rnorm(100)))
#'
#' # Multiple B matrices
#' A <- Diagonal(100)
#' B1 <- rsparsematrix(100, 100, 1, symmetric = TRUE)
#' B2 <- rsparsematrix(100, 100, 1, symmetric = TRUE)
#' W <- WoodburyMatrix(A, B = list(B1, B2))
#' str(solve(W, rnorm(100)))
#' @seealso \linkS4class{WoodburyMatrix}, \link{solve}, \link{instantiate}
#' @references More information on the underlying linear algebra can be found
#' in Harville, D. A. (1997) <doi:10.1007/b98818>.
#' @export
WoodburyMatrix <- function(
  A,
  B,
  U,
  V,
  X,
  O,
  symmetric
) {
  if (!missing(X) && (!missing(U) || !missing(V))) {
    stop('Cannot provide X at the same time as U or V')
  }

  B_parts <- if (!is.list(B)) list(B) else B

  if (missing(symmetric)) {
    if (!missing(U) || !missing(V)) {
      symmetric <- FALSE
    } else if (!missing(X)) {
      symmetric <- TRUE
    } else {
      symmetric <- isSymmetric(A) && all(sapply(B_parts, isSymmetric))
    }
  }

  # Convert everything to lists
  to_list <- function(x, is_missing) {
    if (is_missing) list(Diagonal(nrow(A))) else if (!is.list(x)) list(x) else x
  }
  U_parts <- to_list(U, missing(U))
  V_parts <- to_list(V, missing(V))
  X_parts <- to_list(X, missing(X))

  # Convert to appropriate Matrix subtypes (done before recycling to avoid
  # unnecessarily repeating conversions)
  to_right_Matrix <- function(x) {
    if (is(x, 'diagonalMatrix')) {
      x
    } else {
      as(x, if (symmetric) 'symmetricMatrix' else 'Matrix')
    }
  }
  A <- to_right_Matrix(A)
  B_parts <- lapply(B_parts, to_right_Matrix)
  U_parts <- lapply(U_parts, as, 'Matrix')
  V_parts <- lapply(V_parts, as, 'Matrix')
  X_parts <- lapply(X_parts, as, 'Matrix')

  # Recycle arguments (copy-on-write means arguments aren't actually copied)
  n_parts <- max(sapply(list(B_parts, U_parts, V_parts, X_parts), length))
  B_parts <- .recycle_vector_to(B_parts, n_parts)
  U_parts <- .recycle_vector_to(U_parts, n_parts)
  V_parts <- .recycle_vector_to(V_parts, n_parts)
  X_parts <- .recycle_vector_to(X_parts, n_parts)

  # Check compatibility
  check_dims_match <- function(l, l_dim, r, r_dim) {
    stopifnot(all(sapply(seq_along(l), function(i) {
      l_dim(l[[i]]) == r_dim(r[[i]])
    })))
  }
  check_matches_A <- function(x, x_dim) {
    stopifnot(all(sapply(x, function(U) nrow(A) == x_dim(U))))
  }
  if (!symmetric) {
    check_matches_A(U_parts, nrow)
    check_matches_A(V_parts, ncol)
    check_dims_match(U_parts, ncol, B_parts, nrow)
    check_dims_match(B_parts, ncol, V_parts, nrow)
  } else {
    check_matches_A(X_parts, nrow)
    check_dims_match(X_parts, ncol, B_parts, nrow)
  }

  # Combine parts
  B <- if (n_parts == 1) {
    # NOTE(mgnb): avoid call to .bdiag to avoid potential coercion to sparse
    B_parts[[1]]
  } else {
    .bdiag(B_parts[((seq_len(n_parts) - 1) %% length(B_parts)) + 1])
  }
  U <- if (symmetric) NULL else do.call(cbind, U_parts)
  V <- if (symmetric) NULL else do.call(rbind, V_parts)
  X <- if (!symmetric) NULL else do.call(cbind, X_parts)

  # Object construction
  if (missing(O)) {
    if (!symmetric) {
      O <- B + V %*% A %*% U
    } else {
      rhs <- crossprod(X, A %*% X)
      O <- B + rhs
      if (is(B, 'denseMatrix') || is(rhs, 'denseMatrix')) {
        # HACK(mgnb): Matrix package doesn't coerce to dense matrix in this
        # situation, but it's desirable
        O <- as(O, 'dsyMatrix')
      }
    }
  }
  O <- to_right_Matrix(O)

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
