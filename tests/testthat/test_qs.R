context("Tests that qs function works.")

fit <- qs(hp ~ mpg, data = mtcars)

de <- distributional_effects(fit)
de_mat <- distributional_effects(fit, newdata = head(mtcars))

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
})

testthat::test_that("Prediction functions work", {
  testthat::expect_type(predict(fit), "double")
  testthat::expect_true(is.matrix(predict(fit)))
  testthat::expect_true(is.matrix(predict(fit, newdata = head(mtcars))))
})
