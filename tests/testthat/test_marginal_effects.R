
fit1 <- qs(mpg ~ cyl, data = mtcars, calc_se = F)
fit2 <- qs(mpg ~ cyl + hp, data = mtcars, calc_se = F)
fit3 <- qs(mpg ~ factor(cyl), data = mtcars, calc_se = F)
fit4 <- qs(peri ~ ., data = lapply(rock, function(x) as.numeric(scale(x))),
           algorithm = "br", calc_se = F)

testthat::test_that("Marginal Effects Run OK",{
  testthat::expect_equal(length(marginal_effects(fit1)), 1)
  testthat::expect_equal(length(marginal_effects(fit2)), 2)
  # testthat::expect_equal(attr(marginal_effects(fit1), "type"), "ame")
  testthat::expect_equal(nrow(marginal_effects(fit1)[[1]]), 1)
  testthat::expect_gt(nrow(marginal_effects(fit1, type = "varying")[[1]]), 1)
  testthat::expect_equal(length(marginal_effects(fit4)), ncol(rock) - 1)
  testthat::expect_equal(nrow(marginal_effects(fit3, type = "varying")[[1]]), 2)
})
