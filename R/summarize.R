
#' Generates summary output for quantile and OLS fits
#' @param quant_fit Named rows of quantile beta estimates.
#' @param ols_fit Named rows of OLS beta estimates.
#' @param se Vector containing covariance matrices for both quantile and OLS fits.
#' @param alpha Column vector representing quantiles being estimated.
#' @param varnames vector of variable names for output
#' @param ols_r_squared R-squared for the OLS fit.
#' @return Summary output matrix
summary_output = function(
  quant_fit,
  ols_fit,
  se,
  alpha,
  varnames,
  ols_r_squared) {

  quant_betas <- quant_fit$coef
  pseudo_r <- quant_fit$pseudo_r

  quant_se <- diag(se$quant_cov)^(0.5)
  quant_out_vec <- c()
  quant_out_vec[seq(1,length(quant_betas)*2,2)] <- unlist(quant_betas)
  quant_out_vec[seq(2,length(quant_betas)*2,2)] <- unlist(quant_se)

  num_betas <- length(ols_fit)
  quant_out <- matrix(quant_out_vec, (num_betas*2), length(alpha))
  quant_out <- rbind(quant_out, pseudo_r)

  ols_se <- diag(se$ols_cov)^(0.5)
  ols_out <- c()
  ols_out[seq(1,length(ols_fit)*2,2)] <- unlist(ols_fit)
  ols_out[seq(2,length(ols_fit)*2,2)] <- unlist(ols_se)
  ols_out <- append(ols_out, ols_r_squared)

  final_output <- cbind(quant_out, ols_out)
  colnames(final_output) <- c(alpha, "OLS")

  r_names <- c()
  r_names[seq(1,num_betas*2,2)] <- varnames
  r_names[seq(2,num_betas*2,2)] <- paste(varnames, 'SE', sep = '_')
  r_names <- append(r_names, 'Pseudo-R^2 (Quantile); R^2 (OLS)')
  rownames(final_output) = r_names

  return(final_output)
}
