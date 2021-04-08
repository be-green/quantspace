context("Test that printing helper functions work")

set.seed(10)

fit <- qs(mpg ~ cyl, data = mtcars, calc_se = F)

# The solution should be non-unique
# don't want it to throw a warning
suppressWarnings({
  fit_lasso <- qs(mpg ~ ., data = mtcars, calc_se = F,
                  algorithm = "lasso", control = qs_control(lambda = 1))
})

testthat::test_that("qs print and rounding functions work", {
  testthat::expect_equal(round_if(data.frame(a = "a", n = 2.1), d = 0),
                         data.frame(a = "a", n = 2))
  testthat::expect_output(print(fit))
  testthat::expect_output(print(summary(fit)))
  testthat::expect_output(print(fit_lasso), regexp = "Inter", fixed = T)
  testthat::expect_output(print(fit_lasso), regexp = "^(?!.*(gear))", perl = T)
})
