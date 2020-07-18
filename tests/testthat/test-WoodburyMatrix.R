context('WoodburyMatrix')

load('matrices.rda')

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

  expect_equal(max(abs(instantiate(W) - W_direct)), 0)
  expect_equal(W %*% b, W_direct %*% b)
  expect_equal(determinant(W), determinant(W_direct))
  expect_equal(
    determinant(W, logarithm = FALSE),
    determinant(W_direct, logarithm = FALSE)
  )
  expect_equal(max(abs(solve(W) - solve(W_direct))), 0)
  expect_equal(solve(W, b), solve(W_direct, b))
  expect_equal(max(abs(instantiate(t(W)) - t(W_direct))), 0)
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
}

test_that('GWoodbury with diagonal matrices', {
  n <- 100
  compare_to_direct(Diagonal(n), Diagonal(n, 2))
})

test_that('GWoodbury with general matrices', {
  n <- 100
  compare_to_direct(G_100, G_100)
})

test_that('GWoodbury with general matrices and U and V', {
  n <- 100
  p <- 50
  compare_to_direct(
    G_100,
    G_50,
    U = G_100_50,
    V = t(G_100_50)
  )
})

test_that('SWoodbury operations with negative semidefinite matrices', {
  n <- 100
  compare_to_direct(
    S_100_nd,
    S_100_nd,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with negative semidefinite matrices and X', {
  n <- 100
  p <- 50
  compare_to_direct(
    S_100_nd,
    S_50_nd,
    X = G_100_50,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with positive semidefinite matrices', {
  n <- 100
  compare_to_direct(
    S_100_pd,
    S_100_pd,
    symmetric = TRUE
  )
})

test_that('SWoodbury operations with positive semidefinite matrices and X', {
  n <- 100
  p <- 50
  compare_to_direct(
    S_100_pd,
    S_50_pd,
    X = G_100_50,
    symmetric = TRUE
  )
})

test_that('Creation fails when both U/V and X are provided', {
  n <- 100
  A <- Diagonal(n)
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
