context("Tests that qs function works.")

set.seed(100)
setCores(1)

x = matrix(rnorm(1000), ncol = 2)
y = 1 + 2 * x[,1] - 0.4 * x[,2] + rnorm(nrow(x)) * ( x[,1]) * 4 + rnorm(nrow(x)) * ( x[,2]) *3

test_data = data.frame(y = y, x)

fit <- qs(y ~ X1 + X2, data = head(test_data, 900), parallel = F)

fit_lasso_no_penalty <- qs(y ~ X1 + X2, data = head(test_data, 900),
          parallel = F, scale_x = F,
          algorithm = "lasso",
          method = "br",
          lambda = 0)

fit_br <- qs(y ~ X1 + X2, data = head(test_data, 900),
                           parallel = F,
                           algorithm = "rq.fit.br")


lasso_diff = max(abs(coef(fit_br) - coef(fit_lasso_no_penalty)))


de <- distributional_effects(fit)
de_mat <- distributional_effects(fit, newdata = tail(test_data, 5))

testthat::test_that("S3 Classes inherit properly", {
  testthat::expect_s3_class(fit, "qs")
  testthat::expect_s3_class(summary(fit), "qs_summary")
  testthat::expect_s3_class(plot(de), "ggplot")
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



testthat::test_that("Prediction functions work", {
  testthat::expect_type(predict(fit), "double")
  testthat::expect_true(is.matrix(predict(fit)))
  testthat::expect_true(is.matrix(predict(fit, newdata = tail(test_data, 100))))
})

testthat::test_that("Lasso matches br when not penalized", {
  testthat::expect_equivalent(0, lasso_diff)
})
