if (requireNamespace('lintr', quietly = TRUE)) {
  context('lints')

  test_that('Package passes lintr', {
    lintr::expect_lint_free()
  })
}
