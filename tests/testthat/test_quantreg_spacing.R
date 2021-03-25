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
