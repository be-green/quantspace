#' Compute quantile regressions via quantile spacings
#' @param formula an object of class "formula" (or one that can
#' be coerced to that class): a symbolic description of the model
#' to be fitted.
#' @param data an optional data frame, list or environment (or object
#' coercible by as.data.frame to a data frame) containing the variables in the model.
#' If not found in data, the variables are taken from environment(formula)
#' @param quantiles vector of quantiles to be estimated
#' @param baseline_quantile baseline quantile to measure spacings from (defaults to 0.5)
#' @param se_method one of "delta" or "boot". Delta uses the delta method approximation,
#' while boot bootstraps the standard errors.
#' @importFrom assertthat assert_that
#' @importFrom data.table dcast
#' @export
# qs <- function(formula, data,
#                quantiles = c(0.9, 0.75, 0.5, 0.25, 0.1),
#                baseline_quantile = 0.5,
#                se_method = "delta",
#                ...) {
#
#   assertthat::assert_that(length(baseline_quantile) == 1)
#
#   m <- model.matrix(formula, data)
#   y <- model.response(model.frame(formula, data), type = "numeric")
#
#   quantiles <- unique(quantiles)
#
#   if(!baseline_quantile %in% quantiles){
#     quantiles <- c(quantiles, baseline_quantile)
#   }
#
#   quantiles <- sort(quantiles)
#
#   output <- quantRegSpacing(y,
#                             data = m,
#                             var_names = colnames(m),
#                             alpha = quantiles,
#                             jstar = which(quantiles == baseline_quantile),
#                             ...)
#
#   output$m <- m
#   output$y <- y
#
#   # output <- lapply(output, as.vector)
#
#   # output$se <- ses
#
#   structure(output, class = "qs")
# }

qs <- function(formula, data,
               quantiles = c(0.9, 0.75, 0.5, 0.25, 0.1),
               baseline_quantile = 0.5,
               se_method = "delta",
               trunc = T, small = 1e-3,   weight_vec = NULL,
               outputQuantiles = FALSE, calculateAvgME = FALSE,
               subsamplePct = 0.2, cluster_indices = NULL,
               stratum_indices = NULL, draw_weights = FALSE,
               num_bs = 100, parallelize = FALSE,
               num_cores = 1, save_rows = FALSE, seed = NULL,
               outToCsv = FALSE, csv_filename=NULL,
               ...) {

  assertthat::assert_that(length(baseline_quantile) == 1)

  m <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data), type = "numeric")

  depCol <- y

  quantiles <- unique(quantiles)

  if(!baseline_quantile %in% quantiles){
    quantiles <- c(quantiles, baseline_quantile)
  }

  reg_spec <- denseMatrixToSparse(m)
  reg_spec_var_names <- colnames(m)

  outputQuantiles = T
  calculateAvgME = T

  alpha <- sort(quantiles)

  jstar <- which(quantiles == baseline_quantile)


  message("Calculating initial quantile fit")
  quantreg_fit <- quantRegSpacing(
    dep_col = depCol,
    data = reg_spec,
    var_names = reg_spec_var_names,
    alpha = alpha,
    jstar = jstar,
    outputQuantiles = outputQuantiles,
    calculateAvgME = calculateAvgME,
    ...
  )
  message("Initial quantile fit complete")

  message("Calculating initial OLS fit")
  ptm <- proc.time()
  ols_fit <- ols_sparse_fit(
    a = reg_spec,
    y = depCol)
  message(proc.time() - ptm)
  message("Initial OLS fit complete")

  ols_r_squared <- get_ols_r_squared(
    a = reg_spec,
    y = depCol,
    betas = ols_fit)

  if(subsamplePct > 1) subsamplePct = subsamplePct * 0.01

  message("Starting subsampling")
  se = subsampleStandardErrors(
    dep_col = depCol,
    data = reg_spec,
    var_names = reg_spec_var_names,
    alpha = alpha,
    jstar = jstar,
    M = subsamplePct,
    cluster_indices = cluster_indices,
    stratum_indices = stratum_indices,
    draw_weights = draw_weights,
    num_bs = num_bs,
    parallelize = parallelize,
    num_cores = num_cores,
    trunc = trunc,
    start_model = quantreg_fit$coef,
    small = small,
    save_rows = save_rows)
  message("Subsampling complete")

  rv = list('quantreg_fit' = quantreg_fit,
            'ols_fit' = ols_fit,
            'ols_r_squared' = ols_r_squared,
            'se' = se,
            'specs' = list('formula' = formula,
                           'X' = data,
                           'Y' = depCol,
                           'alpha' = alpha,
                           'jstar' = jstar,
                           'trunc' = trunc,
                           'small' = small,
                           'subsamplePct' = subsamplePct,
                           'cluster_indices'= cluster_indices,
                           'stratum_indices' = stratum_indices,
                           'draw_weights' = draw_weights,
                           'num_bs' = num_bs,
                           'parallelize' = parallelize,
                           'num_cores' = num_cores,
                           'outToCsv' = outToCsv,
                           'csv_filename' = csv_filename))



  structure(rv, class = "qs")
}


summary.qs = function(x, ...){
  # creates a table of summary output for a qs object
  #
  # Args:
  #   qs.class: an object of class qs
  #   createCsv: boolean value, indicating whether csv's should be saved
  #   csvFilename: option value, beginning of name of the csv files
  #
  # Outputs:
  #   matrix of summary statistics for object of class qs
  #   also saves csv files if user turns on this option
  #
  quant_betas <- x$quantreg_fit$coef
  pseudo_r <- x$quantreg_fit$pseudo_r

  quant_se <- diag(x$se$quant_cov)^(0.5)
  quant_out_vec <- c()
  quant_out_vec[seq(1,length(quant_betas)*2,2)] <- unlist(quant_betas)
  quant_out_vec[seq(2,length(quant_betas)*2,2)] <- unlist(quant_se)

  num_betas <- length(x$ols_fit)
  quant_out <- matrix(quant_out_vec, (num_betas*2), length(x$specs$alpha))
  quant_out <- rbind(quant_out, pseudo_r)

  ols_se <- diag(x$se$ols_cov)^(0.5)
  ols_out <- c()
  ols_out[seq(1,length(x$ols_fit)*2,2)] <- unlist(x$ols_fit)
  ols_out[seq(2,length(x$ols_fit)*2,2)] <- unlist(ols_se)
  ols_out <- append(ols_out, x$ols_r_squared)

  final_output <- cbind(quant_out, ols_out)
  colnames(final_output) <- c(alpha, "OLS")

  b_names <- gsub("[0-9\\.]+_",
                  replacement = "",
                  names(
                    head(unlist(x$quantreg_fit$coef),
                         num_betas)
                    )
                  )

  r_names <- c()
  r_names[seq(1,num_betas*2,2)] <- b_names
  r_names[seq(2,num_betas*2,2)] <- paste(b_names, 'se', sep = '_')
  r_names <- append(r_names, 'Pseudo-R^2 (Quantile); R^2 (OLS)')
  rownames(final_output) = r_names
  print(final_output, ...)
}


