
x = matrix(rnorm(1000), ncol = 10)
x <- cbind(1, x, x[,1] + x[,2])


testthat::test_that("Singular design matrix throws an error", {
  testhat::expect_error(ensureSpecFullRank(x))

})

