context("Tests that qs & distributional_effects work.")

set.seed(100)
setCores(1)

x = matrix(rnorm(1000), ncol = 2)
y = 1 + 2 * x[,1] - 0.4 * x[,2] + rnorm(nrow(x)) * ( x[,1]) * 4 + rnorm(nrow(x)) * ( x[,2]) *3

test_data = data.frame(y = y, x)

fit <- qs(y ~ X1 + X2, data = head(test_data, 900), algorithm = "agd", parallel = F)

de <- distributional_effects(fit)
de_mat <- distributional_effects(fit, newdata = tail(test_data, 5))

fit_no_se <- qs(mpg ~ cyl, data = mtcars, parallel = F, calc_se = F)
fit_with_me <- qs(mpg ~ cyl, data = mtcars, parallel = F, calc_se = F,
                  control = qs_control(calc_avg_me = T))

testthat::test_that("S3 Classes inherit properly", {
  testthat::expect_s3_class(fit, "qs")
  testthat::expect_s3_class(summary(fit), "qs_summary")
  testthat::expect_s3_class(plot(de), "ggplot")
})

testthat::test_that("qs errors and warnings work", {
  testthat::expect_error(qs(y ~ X1, data = test_data, algorithm = "UNIMPLEMENTED_ALGORITHM"))
  testthat::expect_warning(qs(y ~ X1, data = head(test_data, 3),
                              std_err_control = se_control(se_method = "bootstrap")))
  testthat::expect_error(qs(y ~ X1, data = head(test_data, 1), parallel = F))
  testthat::expect_error(qs(y ~ X1, data = head(test_data, 1), parallel = F,
                            std_err_control = se_control(subsample_percent = 10)))
  testthat::expect_error(qs(y ~ X1, data = head(test_data, 1), parallel = F,
                            std_err_control = se_control(subsample_percent = -10)))
  testthat::expect_error(qs(y ~ X1, data = head(test_data, 1), parallel = F,
                            se_method = "UNIMPLEMENTED_ALGORITHM"))
})

testthat::test_that("calc_se option doesn't calculate ses when set to false",{
  testthat::expect_true(all(is.na(fit_no_se$se$quant_cov)))
})

test_baseline = qs(mpg ~ cyl, data = mtcars, calc_se = F, baseline_quantile = 0.5, quantiles = c(0.25, 0.75))

testthat::test_that("Baseline quantile ends up in the mix if it isn't already present in quantiles", {
  testthat::expect_true(0.5 %in% test_baseline$specs$alpha)
})


testthat::test_that("qs returns marginal effects if required", {
  testthat::expect_true(is.numeric(fit_with_me$quantreg_fit$me))
})

testthat::test_that("Distributional effect functions work", {
  testthat::expect_type(de$r(10), "double")
  testthat::expect_length(de$r(10), 10)
  testthat::expect_type(de$q(0.8), "double")
  testthat::expect_error(de$q(10))
  testthat::expect_type(de_mat$r(10), "double")
  testthat::expect_equal(ncol(de_mat$r(10)), 10)
  testthat::expect_type(de_mat$q(0.8), "double")
  testthat::expect_error(de_mat$q(10))
  testthat::expect_true(inherits(plot(de), "gg"))
  testthat::expect_true(inherits(plot(de_mat), "gg"))
})

fit_no_quantiles <- qs(y ~ X1, data = head(test_data, 100), output_quantiles = F, calc_se = F)

testthat::test_that("Prediction functions work", {
  testthat::expect_type(predict(fit), "double")
  testthat::expect_true(is.matrix(predict(fit)))
  testthat::expect_true(is.matrix(predict(fit, newdata = tail(test_data, 100))))
  testthat::expect_true(is.matrix(predict(fit_no_quantiles)))
  testthat::expect_true(is.matrix(predict(fit_no_quantiles, newdata = tail(test_data, 100))))
})

testthat::test_that("NA if null works",{
  testthat::expect_equal(na_if_null(NA), NA)
  testthat::expect_equal(na_if_null(NULL), NA)
})
