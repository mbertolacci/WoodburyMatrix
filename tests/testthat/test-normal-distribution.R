context('normal distribution')

test_that('dwnorm returns the correct value for known examples', {
  n <- 100
  W <- WoodburyMatrix(Diagonal(n), Diagonal(n))

  true_with_0 <- -0.5 * (n * log(2 * pi) + n * log(2))
  true_with_1 <- true_with_0 - 0.5 * n / 2

  expect_equal(
    dwnorm(rep(1, n), covariance = W, log = TRUE),
    true_with_1
  )
  expect_equal(
    dwnorm(rep(1, n), covariance = W, log = FALSE),
    exp(true_with_1)
  )
  expect_equal(
    dwnorm(rbind(rep(0, n), rep(1, n)), covariance = W, log = TRUE),
    c(true_with_0, true_with_1)
  )
})

test_that('rwnorm return value has correct dimensions', {
  n <- 100

  check_shapes <- function(W) {
    expect_length(rwnorm(1, covariance = W), n)
    expect_length(rwnorm(1, 1, W), n)
    expect_equal(dim(rwnorm(2, 1, W)), c(2, n))
  }

  D <- Diagonal(n)
  check_shapes(WoodburyMatrix(D, D))
  D_dense <- as(as(D, 'dgeMatrix'), 'dsyMatrix')
  check_shapes(WoodburyMatrix(D_dense, D_dense))
})
