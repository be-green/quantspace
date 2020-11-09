context("Tests that function works.")

fit <- qs(hp ~ mpg, data = mtcars, quantiles = c(0.3, 0.5, 0.7))

testthat::test_that("Regression runs", {
  testthat::expect_equal(class(fit), "qs")
})
