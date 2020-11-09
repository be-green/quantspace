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
#' @param trunc whether to truncate small values
#' @param small level of "small" values to guarentee numerical stability
#' @param weight_vec vector of weights for weighted quantile regression
#' @param outputQuantiles whether to output the fitted quantiles
#' @param calculateAvgME whether to calculated average marginal effects
#' @param subsamplePct percent to subsample for standard error calculations
#' @param cluster_indices index of clusters for clustered standard errors
#' @param stratum_indices index of strata for clustered standard errors
#' @param draw_weights whether to use random exponential weights for bootstrap
#' @param num_bs number of bootstrap draws
#' @param parallel whether to run bootstrap in parallel
#' @param num_cores number of cores to use (defaults to setting from `getOption(mc.cores)`)
#' @param save_rows whether to save rows, see [subsampleStandardErrors]
#' @param seed what seed to use for replicable RNG
#' @param outToCsv whether to save fit to a csv file
#' @param csv_filename filename of csv outfile
#' @param ... additional arguments, ignored for now
#' @importFrom assertthat assert_that
#' @importFrom data.table dcast
#' @importFrom SparseM model.matrix
#' @importFrom stats model.frame
#' @importFrom SparseM model.response
#' @export
qs <- function(formula, data,
               quantiles = c(0.9, 0.75, 0.5, 0.25, 0.1),
               baseline_quantile = 0.5,
               se_method = "boot",
               trunc = T, small = 1e-3, weight_vec = NULL,
               outputQuantiles = TRUE, calculateAvgME = TRUE,
               subsamplePct = 0.2, cluster_indices = NULL,
               stratum_indices = NULL, draw_weights = FALSE,
               num_bs = 100, parallel = TRUE,
               num_cores = NULL, save_rows = FALSE, seed = NULL,
               outToCsv = FALSE, csv_filename = NULL,
               ...) {



  assertthat::assert_that(length(baseline_quantile) == 1)

  if(is.null(num_cores)) {
    num_cores <- getCores()
  }

  m <- SparseM::model.matrix(formula, data)
  y <- SparseM::model.response(stats::model.frame(formula, data), type = "numeric")

  depCol <- y

  quantiles <- unique(quantiles)

  if(!baseline_quantile %in% quantiles){
    quantiles <- c(quantiles, baseline_quantile)
  }

  reg_spec <- denseMatrixToSparse(m)
  reg_spec_var_names <- colnames(m)

  alpha <- sort(quantiles)

  jstar <- which(quantiles == baseline_quantile)

  # message("Calculating initial quantile fit")
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

  ols_fit <- ols_sparse_fit(
    a = reg_spec,
    y = depCol)

  ols_r_squared <- get_ols_r_squared(
    a = reg_spec,
    y = depCol,
    betas = ols_fit)

  if(subsamplePct > 1) subsamplePct = subsamplePct * 0.01

  if(se_method == "boot") {
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
      parallel = parallel,
      num_cores = num_cores,
      trunc = trunc,
      start_model = quantreg_fit$coef,
      small = small,
      save_rows = save_rows)
  } else {
    stop("This method is not yet implemented or integrated with qs.")
  }


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
                           'parallel' = parallel,
                           'num_cores' = num_cores,
                           'outToCsv' = outToCsv,
                           'csv_filename' = csv_filename))

  structure(rv, class = "qs")
}

#' creates a table of summary output for a qs object
#' @param x an object of class qs
#' @param ... other arguments to pass to summary
#' @importFrom utils head
summary.qs = function(x, ...){

  quant_betas <- x$quantreg_fit$coef
  pseudo_r <- x$quantreg_fit$pseudo_r

  quant_se <- diag(x$se$quant_cov)^(0.5)
  quant_out_vec <- c()
  quant_out_vec[seq(1,length(quant_betas)*2,2)] <- unlist(quant_betas)
  quant_out_vec[seq(2,length(quant_betas)*2,2)] <- unlist(quant_se)

  num_betas <- length(x$ols_fit)
  quant_out <- matrix(quant_out_vec, (num_betas*2), length(x$specs$alpha))

  ols_se <- sqrt(diag(x$se$ols_cov))
  ols_out <- rbind(x$ols_fit, ols_se)

  coef_mat <- cbind(quant_out)
  colnames(coef_mat) <- c(x$specs$alpha)

  b_names <- gsub("[0-9\\.]+_",
                  replacement = "",
                  names(
                    utils::head(unlist(x$quantreg_fit$coef),
                         num_betas)
                    )
                  )


  r_names <- c()
  r_names[seq(1,num_betas*2,2)] <- b_names
  r_names[seq(2,num_betas*2,2)] <- paste(b_names, 'se', sep = '_')
  rownames(coef_mat) <- r_names

  final_output <- list(
    q_coefs = coef_mat,
    ols_coefs = ols_out,
    R2 = list(psuedo_r = pseudo_r,
              ols_r_squared = x$ols_r_squared)
  )

  structure(final_output, class = "qs_summary")
}
#' Print qs summary
#' @param x fit from qs
#' @param d number of digits to print
#' @param ... additional arguments, ignored for now
#' @export
print.qs <- function(x, d = 4, ...) {

  # really poor solution, going with it for now
  x <- summary(x, ...)

  q_coefs <- x$q_coefs
  nms <- rownames(q_coefs)
  nms[which(1:length(nms) %% 2 == 0)] <- ""
  rownames(q_coefs) <- nms
  q_coefs <- make_se_mat(q_coefs, d)

  cat("Quantile Coefficients:\n",
      paste0(q_coefs, "\n\n"))

  o_coefs <- x$ols_coefs
  rownames(o_coefs) <- rep("", nrow(o_coefs))
  o_coefs <- make_se_mat(o_coefs, d)

  cat(
    "OLS Coefficients:\n",
    o_coefs)

}

#' Make matrix into a "standard error" matrix
#' @param mat matrix to pass
#' @param d number of significant digits
#' @importFrom testthat capture_output
#' @importFrom stringr str_replace_all
make_se_mat <- function(mat, d) {
  mat <- signif(mat, d)

  msg <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for(i in 1:nrow(mat)) {
    if(i %% 2 == 0) {
      this_row <- paste0(
        sapply(mat[i,], function(x) paste0("(",x,")"))
      )
    } else {
      this_row <- sapply(mat[i,], as.character)
    }

    msg[i,] <- this_row

  }
  colnames(msg) <- colnames(mat)
  rownames(msg) <- rownames(mat)
  msg <- testthat::capture_output(print(msg))
  stringr::str_replace_all(msg, "\"", " ")
}

#' Print qs summary
#' @param x fit from qs
#' @param d number of digits to print
#' @param ... additional arguments, ignored for now
#' @export
print.qs_summary <- function(x, d = 4, ...) {
  args <- list(...)

  psuedo_r2 <- signif(x$R2$psuedo_r, d)
  psuedo_message <- paste0(paste0("    ",
                                      colnames(x$q_coefs),
                                      ":\t",
                                      psuedo_r2), collapse = "\n")

  o_coefs <- x$ols_coefs
  rownames(o_coefs) <- rep("", nrow(o_coefs))
  o_coefs <- make_se_mat(o_coefs, d)

  q_coefs <- x$q_coefs
  nms <- rownames(q_coefs)
  nms[which(1:length(nms) %% 2 == 0)] <- ""
  rownames(q_coefs) <- nms
  q_coefs <- make_se_mat(q_coefs, d)

  cat("Quantile Coefficients:\n",
      paste0(q_coefs, "\n\n"))

  cat(paste0("Quantile Psuedo-R2:\n",
       psuedo_message, "\n\n"))
  cat(
    "OLS Coefficients:\n",
    o_coefs, "\n")
  cat("R2: ", signif(x$R2$ols_r_squared, d))
}

