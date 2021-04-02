context("Test Quantile Regression Algorithm Interface")

X <- SparseM::model.matrix(mpg ~ cyl, data = mtcars)
y <- SparseM::model.response(stats::model.frame(mpg ~ cyl, data = mtcars), type = "numeric")

reg_spec <- denseMatrixToSparse(X)

star_model = rq.fit.sfn_start_val(
  X = reg_spec,
  y = y,
  tau = 0.5,
  weight_vec = NULL)

test_star_model <- fitQuantileRegression(
  X = reg_spec,
  y = y,
  tau = 0.5,
  weight_vec = NULL,
  algorithm = "rq.fit.sfn_start_val"
)

test_that("Generic Algorithm Interface Matches Specific Output", {
  testthat::expect_equivalent(star_model, test_star_model)
  testthat::expect_equal(class(X), c("matrix", "array"))
  testthat::expect_true(inherits(reg_spec, "matrix.csr"))
})


x = matrix(rnorm(1000), ncol = 2)
y = 1 + 2 * x[,1] - 0.4 * x[,2] + rnorm(nrow(x)) * ( x[,1]) * 4 + rnorm(nrow(x)) * ( x[,2]) *3

test_data = data.frame(y = y, x)

fit_lasso_no_penalty <- qs(y ~ X1 + X2, data = head(test_data, 900),
                           parallel = F, scale_x = F,
                           algorithm = "lasso",
                           method = "br",
                           lambda = 0)

fit_br <- qs(y ~ X1 + X2, data = head(test_data, 900),
             parallel = F,
             algorithm = "rq.fit.br")


lasso_diff = max(abs(coef(fit_br) - coef(fit_lasso_no_penalty)))

testthat::test_that("Lasso matches br when not penalized", {
  testthat::expect_equivalent(0, lasso_diff)
})

