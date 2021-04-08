context("Test utility for printing quantile regression warnings")

testthat::test_that("error printing utility works for quantreg_spacing", {
  testthat::expect_warning(printWarnings(list(ierr="test", it = 100)))
  testthat::expect_warning(printWarnings(list(ierr="test", it = 2)))
})
