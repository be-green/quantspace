context("Test standard errors")

suppressWarnings({
  boot_fit <- qs(mpg ~ hp, data = mtcars, se_method = "bootstrap")
  wboot_fit <- qs(mpg ~ hp, data = mtcars, se_method = "weighted_bootstrap")
  subsample_fit <- qs(mpg ~ hp, data = mtcars, se_method = "subsample")
  custom_fit <- qs(mpg ~ hp, data = mtcars, se_method = "resample_qs",
                   draw_weights = T, sampling_method = "subsample",
                   subsample_percent = 0.3)
})


test_that("Standard Error Algorithms Run",{
  testthat::expect_true(inherits(boot_fit, "qs"))
  testthat::expect_true(inherits(wboot_fit, "qs"))
  testthat::expect_true(inherits(subsample_fit, "qs"))
  testthat::expect_true(inherits(custom_fit, "qs"))
})
