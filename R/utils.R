.try_Cholesky <- function(...) {
  tryCatch(
    Cholesky(...),
    error = function(e) {},
    warning = function(w) {}
  )
}

.chol_Linv_P <- function(a, b) {
  solve(a, solve(a, b, system = 'P'), system = 'L')
}
