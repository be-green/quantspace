
fit1 <- qs(mpg ~ cyl, data = mtcars, calc_se = F)
fit2 <- qs(mpg ~ cyl + hp, data = mtcars, calc_se = F)
fit3 <- qs(mpg ~ factor(cyl), data = mtcars, calc_se = F)
fit4 <- qs(peri ~ ., data = lapply(rock, function(x) as.numeric(scale(x))),
           algorithm = "br", calc_se = F)
suppressWarnings({

fit5 <- qs(peri ~ ., data = lapply(rock, function(x) as.numeric(scale(x))),
           algorithm = "br", calc_se = T, std_err_control = se_control(se_method = "bootstrap", num_bs = 10))

})

testthat::test_that("Marginal Effects Run OK",{
  testthat::expect_equal(nrow(marginal_effects(fit1)$avgME), 1)
  testthat::expect_equal(nrow(marginal_effects(fit2)$avgME), 2)
  # testthat::expect_equal(attr(marginal_effects(fit1), "type"), "ame")
  testthat::expect_equal(nrow(marginal_effects(fit1)[[1]]), 1)
  testthat::expect_equal(nrow(marginal_effects(fit4)$avgME), ncol(rock) - 1)
  testthat::expect_equal(length(which(is.na(marginal_effects(fit5)$avgME_se))), 0)
})
