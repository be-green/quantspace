set.seed(42)

x = matrix(rnorm(1000), ncol = 10)
x <- cbind(1, x, x[,1] + x[,2])
x <- cbind(x, rnorm(nrow(x)))
x <- cbind(x, x[,3] * 0.2 - x[,4] * 1.3)

non_colinear_x = matrix(rnorm(1000), ncol = 10)

test_that("Singular design matrix throws an error", {
  expect_warning(new_x <- ensureSpecFullRank(x, paste0("X",1:ncol(x))))
  expect_equal(ncol(new_x$spec_matrix), 12)
  expect_equal(ncol(ensureSpecFullRank(non_colinear_x,
                                       paste0("X",1:ncol(non_colinear_x)))$spec_matrix),
               10)
})

