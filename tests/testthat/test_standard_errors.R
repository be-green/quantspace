set.seed(100)

test_that( "Singularity errors in bootstrapping become warnings", {
  expect_warning(qs(mpg ~ hp, data = mtcars,
                    std_err_control = se_control(se_method = "bootstrap"),
                    parallel = F))
})

X <- matrix(rt(10000, 1.5), ncol = 10)
y <- 0.2 + X %*% rnorm(ncol(X)) + exp(rnorm(nrow(X)))
# these will all generate warnings, ignore it for now

data = data.frame(y = y, X)

suppressWarnings({
  boot_fit <- qs(y ~ X, data = data,
                 std_err_control = se_control(se_method = "bootstrap"),
                 parallel = T)

  wboot_fit <- qs(y ~ X, data = data,
                  std_err_control = se_control(se_method = "weighted_bootstrap"),
                  parallel = T)

  subsample_fit <- qs(y ~ X, data = data,
                      std_err_control = se_control(se_method = "subsample"),
                      parallel = T)

  custom_fit <- qs(y ~ X, data = data,
                   std_err_control = se_control(se_method = "resample_qs",
                   draw_weights = T, sampling_method = "subsample",
                   subsample_percent = 0.3),
                   parallel = T)

  cluster_fit <- qs(mpg ~ hp, data = mtcars,
                    std_err_control = se_control(se_method = "bootstrap"),
                    cluster_formula = ~ cyl)
})


test_that("Standard Error Algorithms Run",{
  testthat::expect_true(inherits(boot_fit, "qs"))
  testthat::expect_true(inherits(wboot_fit, "qs"))
  testthat::expect_true(inherits(subsample_fit, "qs"))
  testthat::expect_true(inherits(custom_fit, "qs"))
  testthat::expect_true(inherits(cluster_fit, "qs"))
})
