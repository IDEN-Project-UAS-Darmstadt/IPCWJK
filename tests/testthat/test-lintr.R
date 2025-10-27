test_that("Code is lint free", {
  if (requireNamespace("lintr", quietly = TRUE)) {
    lintr::expect_lint_free()
  }
})
