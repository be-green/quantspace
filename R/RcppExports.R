# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Drop Colinear Columns from dense matrix
#' @param X matrix to drop colinear columns from
#' @param tol tolerance for when pivot is zero in rank calculations
#' @export
qr_drop_colinear_columns <- function(X, tol = 0.000000001) {
    .Call(`_quantspace_qr_drop_colinear_columns`, X, tol)
}

#' Drop Colinear Columns from sparse matrix
#' @param X matrix to drop colinear columns from
#' @param tol tolerance for when pivot is zero in rank calculations
#' @export
sparse_qr_drop_colinear_columns <- function(X, tol = 0.000000001) {
    .Call(`_quantspace_sparse_qr_drop_colinear_columns`, X, tol)
}

checkfun <- function(res, tau) {
    .Call(`_quantspace_checkfun`, res, tau)
}

parallelVectorCheckFun <- function(x, tau) {
    .Call(`_quantspace_parallelVectorCheckFun`, x, tau)
}

reorder_columns <- function(X, intercept) {
    invisible(.Call(`_quantspace_reorder_columns`, X, intercept))
}

fast_rexp <- function(n) {
    .Call(`_quantspace_fast_rexp`, n)
}

update_huber_grad <- function(X_t, res, derivs, grad, tau, mu, n, one_over_n) {
    invisible(.Call(`_quantspace_update_huber_grad`, X_t, res, derivs, grad, tau, mu, n, one_over_n))
}

z_score <- function(X, colwise_avg_x, colwise_sd_x, p) {
    invisible(.Call(`_quantspace_z_score`, X, colwise_avg_x, colwise_sd_x, p))
}

arma_qr_drop_colinear_columns <- function(X) {
    .Call(`_quantspace_arma_qr_drop_colinear_columns`, X)
}

huber_grad_descent <- function(y, X, X_t, beta, grad, derivs, tau, n, one_over_n, p, maxiter, mu, beta_tol, check_tol, min_delta = 1e-10) {
    .Call(`_quantspace_huber_grad_descent`, y, X, X_t, beta, grad, derivs, tau, n, one_over_n, p, maxiter, mu, beta_tol, check_tol, min_delta)
}

#' Compute quantile regression via accelerated gradient descent using
#' Huber approximation, warm start based on data subset
#' @param y outcome vector
#' @param X design matrix
#' @param X_sub subset of X matrix to use for "warm start" regression
#' @param y_sub subset of y to use for "warm start" regression
#' @param tau target quantile
#' @param init_beta initial guess at beta
#' @param intercept location of the intercept column, using R's indexing
#' @param num_samples number of samples used for subset of matrix used for warm start
#' @param mu neighborhood over which to smooth
#' @param maxiter maximum number of iterations to run
#' @param check_tol loss function change tolerance for early stopping
#' @param beta_tol tolerance for largest element of gradient, used
#' for early stopping
#' @param warm_start integer indicating whether to "warm up" on a subsample
#' of the data
#' @param scale whether to scale x & y variables
#' @param lambda optional lasso penalty weight
#' @param min_delta smallest allowed step size for gradient descent
#' @export
fit_approx_quantile_model <- function(X, y, X_sub, y_sub, tau, init_beta, mu = 1e-15, maxiter = 100000L, beta_tol = 1e-4, check_tol = 1e-6, intercept = 1L, num_samples = 1000, warm_start = 1L, scale = 1L, lambda = 0, min_delta = 1e-10) {
    .Call(`_quantspace_fit_approx_quantile_model`, X, y, X_sub, y_sub, tau, init_beta, mu, maxiter, beta_tol, check_tol, intercept, num_samples, warm_start, scale, lambda, min_delta)
}

#' Compute quantile regression via accelerated gradient descent using
#' Huber approximation, warm start based on data subset
#' @param y outcome vector
#' @param X design matrix
#' @param X_sub subset of X matrix to use for "warm start" regression
#' @param y_sub subset of y to use for "warm start" regression
#' @param tau target quantile
#' @param init_beta initial guess at beta
#' @param intercept location of the intercept column, using R's indexing
#' @param num_samples number of samples used for subset of matrix used for warm start
#' @param mu neighborhood over which to smooth
#' @param maxiter maximum number of iterations to run
#' @param check_tol loss function change tolerance for early stopping
#' @param beta_tol tolerance for largest element of gradient, used
#' for early stopping
#' @param warm_start integer indicating whether to "warm up" on a subsample
#' of the data
#' @param scale whether to scale x & y variables
#' @export
fit_penalize_approx_quantile_model <- function(X, y, X_sub, y_sub, tau, init_beta, mu = 1e-15, maxiter = 100000L, beta_tol = 1e-4, check_tol = 1e-6, intercept = 1L, num_samples = 1000, warm_start = 1L, scale = 1L) {
    .Call(`_quantspace_fit_penalize_approx_quantile_model`, X, y, X_sub, y_sub, tau, init_beta, mu, maxiter, beta_tol, check_tol, intercept, num_samples, warm_start, scale)
}

#' Glob observations w/ residuals above a certain magnitude
#' @param X design matrix
#' @param r vector of residuals
#' @param thresh threshold to use when computing globs
#' @export
glob_obs_mat <- function(X, r, thresh) {
    .Call(`_quantspace_glob_obs_mat`, X, r, thresh)
}

#' Glob observations w/ residuals above a certain magnitude
#' @param y design matrix
#' @param r vector of residuals
#' @param thresh threshold to use when computing globs
glob_obs_vec <- function(y, r, thresh) {
    .Call(`_quantspace_glob_obs_vec`, y, r, thresh)
}

