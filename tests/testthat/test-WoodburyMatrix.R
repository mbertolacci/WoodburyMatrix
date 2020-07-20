context('WoodburyMatrix')

load('matrices.rda')

expect_Matrix_equal <- function(a, b) {
  expect_equal(max(abs(a - b)), 0)
}

compare_to_direct <- function(
  A,
  B,
  U = Diagonal(nrow(B)),
  V = Diagonal(ncol(B)),
  X = Diagonal(nrow(B)),
  symmetric = FALSE
) {
  if (symmetric) {
    W <- WoodburyMatrix(A, B, X = X, symmetric = TRUE)
    W_direct <- forceSymmetric(solve(A) + X %*% solve(B) %*% t(X))
  } else {
    W <- WoodburyMatrix(A, B, U = U, V = V, symmetric = FALSE)
    W_direct <- solve(A) + U %*% solve(B) %*% V
  }
  b <- rnorm(nrow(A))

  expect_Matrix_equal(instantiate(W), W_direct)
  expect_equal(W %*% b, W_direct %*% b)
  expect_equal(determinant(W), determinant(W_direct))
  expect_equal(
    determinant(W, logarithm = FALSE),
    determinant(W_direct, logarithm = FALSE)
  )
  expect_Matrix_equal(solve(W), solve(W_direct))
  expect_equal(solve(W, b), solve(W_direct, b))
  expect_Matrix_equal(instantiate(t(W)), t(W_direct))
  if (symmetric) {
    expect_equal(
      mahalanobis(b, center = 0, cov = W),
      mahalanobis(b, center = 0, cov = W_direct)
    )
    expect_equal(
      mahalanobis(b, center = 1, cov = W),
      mahalanobis(b, center = 1, cov = W_direct)
    )
  }
  expect_equal(isSymmetric(W), symmetric)
}

test_that('GWoodbury with diagonal matrices', {
  n <- 100
  compare_to_direct(Diagonal(n), Diagonal(n, 2))
})

test_that('GWoodbury with general matrices', {
  compare_to_direct(G_100, G_100)
})

test_that('GWoodbury with general matrices and U and V', {
  compare_to_direct(
    G_100,
    G_50,
    U = G_100_50,
    V = t(G_100_50)
  )
})

test_that('SWoodbury operations with negative semidefinite matrices', {
  compare_to_direct(
    S_100_nd,
    S_100_nd,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with negative semidefinite matrices and X', {
  compare_to_direct(
    S_100_nd,
    S_50_nd,
    X = G_100_50,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with positive semidefinite matrices', {
  compare_to_direct(
    S_100_pd,
    S_100_pd,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with positive semidefinite matrices and X', {
  compare_to_direct(
    S_100_pd,
    S_50_pd,
    X = G_100_50,
    symmetric = TRUE
  )
})

test_that('Creation fails when both U/V and X are provided', {
  A <- Diagonal(100)
  expect_error(WoodburyMatrix(A, A, U = A, X = A))
  expect_error(WoodburyMatrix(A, A, V = A, X = A))
  expect_error(WoodburyMatrix(A, A, U = A, V = A, X = A))
})

test_that('Automatic symmetry detection works', {
  n <- 100

  # Matrix type that is definitionally symmetric
  M_ss <- Diagonal(n)
  expect_is(WoodburyMatrix(M_ss, M_ss), 'SWoodburyMatrix')

  # Numerically symmetric
  M_sg <- as(rsparsematrix(n, n, 0.25, symmetric = TRUE), 'dgCMatrix')
  expect_is(WoodburyMatrix(M_sg, M_sg), 'SWoodburyMatrix')

  # Numerically non-symmetric
  M_gg <- rsparsematrix(n, n, 0.25)
  expect_is(WoodburyMatrix(M_ss, M_gg), 'GWoodburyMatrix')
  expect_is(WoodburyMatrix(M_gg, M_gg), 'GWoodburyMatrix')

  # Providing U, V or X forces the issue
  expect_is(WoodburyMatrix(M_ss, M_ss, X = M_ss), 'SWoodburyMatrix')
  expect_is(WoodburyMatrix(M_ss, M_ss, U = M_ss), 'GWoodburyMatrix')
  expect_is(WoodburyMatrix(M_ss, M_ss, V = M_ss), 'GWoodburyMatrix')

  # Should throw an error if non-symmetric matix provided when symmetry = TRUE
  expect_error(WoodburyMatrix(M_gg, M_gg, symmetry = TRUE))
})

test_that('Multiple B matrices can be provided', {
  W1 <- WoodburyMatrix(G_100, list(G_100, 2 * G_100))
  expect_is(W1, 'GWoodburyMatrix')
  expect_Matrix_equal(W1@B, bdiag(G_100, 2 * G_100))
  expect_Matrix_equal(W1@U, cbind(Diagonal(100), Diagonal(100)))
  expect_Matrix_equal(W1@V, rbind(Diagonal(100), Diagonal(100)))

  W2 <- WoodburyMatrix(S_100_pd, list(S_100_pd, 2 * S_100_pd))
  expect_is(W2, 'SWoodburyMatrix')
  expect_Matrix_equal(W2@B, bdiag(S_100_pd, 2 * S_100_pd))
  expect_Matrix_equal(W2@X, cbind(Diagonal(100), Diagonal(100)))
})

test_that('Multiple U matrices can be provided', {
  W <- WoodburyMatrix(G_100, G_100, U = list(Diagonal(100), 2 * Diagonal(100)))
  expect_Matrix_equal(W@U, cbind(Diagonal(100), 2 * Diagonal(100)))
})

test_that('Argument recycling', {
  # Recycling of B works
  W1 <- WoodburyMatrix(G_100, G_100, U = list(Diagonal(100), 2 * Diagonal(100)))
  expect_Matrix_equal(W1@B, bdiag(G_100, G_100))

  # Recycling of U works
  WW <- WoodburyMatrix(G_100, list(G_100, G_100), U = 2 * Diagonal(100))
  expect_Matrix_equal(WW@U, cbind(2 * Diagonal(100), 2 * Diagonal(100)))
})

test_that('O is dense if appropriate', {
  D <- Diagonal(100)
  D_dense <- as(as(D, 'dgeMatrix'), 'dsyMatrix')
  expect_is(WoodburyMatrix(D, D)@O, 'sparseMatrix')
  expect_is(WoodburyMatrix(D_dense, D)@O, 'denseMatrix')
  expect_is(WoodburyMatrix(D, D_dense)@O, 'denseMatrix')
  expect_is(WoodburyMatrix(D_dense, D_dense)@O, 'denseMatrix')
})
