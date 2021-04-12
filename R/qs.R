#' Check if algorithm exists
#' @param algorithm algorithm to check
check_algorithm <- function(algorithm) {
  if(algorithm == "sfn") {
    algorithm = "rq.fit.sfn"
  } else if (algorithm == "lasso") {
    algorithm = "rq.fit.lasso"
  } else if(algorithm == "post_lasso") {
    algorithm = "rq.fit.post_lasso"
  } else if(algorithm == "br") {
    algorithm = "rq.fit.br"
  } else {
    if(!exists(algorithm)) {
      stop(paste0("Algorithm not implemented in quantspace, and not a function" ,
                  " available in the current environment."))
    }
  }
  algorithm
}




#' Control standard_errors parameters
#' @details se_control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`standard_errors`].
#' @param se_method Method to use for standard errors, either "weighted_bootstrap",
#' "subsample", "bootstrap" or "rq_sample" (the last one is fully customized per the other se_control parameters)
#' along with a specified subsampling method and
#' subsample percent. If specifying "custom", must also specify `subsampling_percent` and
#' `draw_weights`. If you specify "subsample", subsampling percent defaults to 0.2, but can be
#' changed.
#' @param num_bs  Number of bootstrap iterations to use, defaults to 100.
#' @param subsample_percent A number between 0 and one, specifying the percent of the data to subsample for standard error calculations
#' @param draw_weights Whether to use random exponential weights for bootstrap, either TRUE or FALSE
#' @param sampling_method One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' @param ... Other arguments, ignored for now
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and bootstrapRows
#' samples with replacement.
#' @param parallel_thresh threshold for when to use parallel, based on `nrow(X) * num_bs`. To always
#' use parallel, set to 0. Designed to avoid slow overhead when problem is relatively small.
#' @export
se_control <- function(se_method = "subsample",
                       subsample_percent = 0.2,
                       num_bs = 100,
                       draw_weights = F,
                       sampling_method = "subsampleRows",
                       parallel_thresh = 100000,
                       ...) {
  list(
    se_method = se_method,
    subsample_percent = subsample_percent,
    num_bs = num_bs,
    draw_weights = draw_weights,
    sampling_method = sampling_method,
    parallel_thresh = parallel_thresh,
    ...)
}

#' Control quantreg_spacing parameters
#' @details control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`qs`]. This is set via the function [`qs_control`],
#' which returns a named list, with elements including:
#' @param trunc whether to truncate residual values below the argument "small"
#' @param small level of "small" values to guarentee numerical stability. If not specified, set dynamically based on the standard deviation of the outcome variable.
#' @param lambda For penalized regression, you can specify a level of lambda which will weight the penalty. If not set, will be determined based on 10-fold cross-validation.
#' @param output_quantiles whether to save fitted quantiles as part of the function output
#' @param calc_avg_me whether to return average marginal effects as part of the fitted object
#' @param ... Other arguments, ignored for now
#' @export
qs_control <- function(trunc = TRUE,
                       small = NULL,
                       lambda = NULL,
                       output_quantiles = TRUE,
                       calc_avg_me = FALSE,
                       ...) {

  list(trunc = trunc,
       small = small,
       lambda = lambda,
       output_quantiles = output_quantiles,
       calc_avg_me = calc_avg_me,
       ...)
}


#' Compute quantile regressions via quantile spacings
#' @param formula an object of class "formula" (or one that can
#' be coerced to that class): a symbolic description of the model
#' to be fitted.
#' @param data an optional data frame, list or environment (or object
#' coercible by as.data.frame to a data frame) containing the variables in the model.
#' If not found in data, the variables are taken from environment(formula)
#' @param quantiles vector of quantiles to be estimated
#' @param baseline_quantile baseline quantile to measure spacings from (defaults to 0.5)
#' @param calc_se boolean, whether or not to calculate standard errors
#' @param weights optional vector of weights for weighted quantile regression
#' @param cluster_formula formula (e.g. ~X1 + X2) giving the clustering formula
#' @param parallel whether to run bootstrap in parallel
#' @param seed what seed to use for replicable RNG
#' @param algorithm What algorithm to use for fitting underlying regressions. Either one of "sfn", "br", "lasso", "post_lasso", or a function name which estimates quantiles. See details.
#' @param control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`qs`]. This is set via the function [`qs_control`],
#' which returns a named list, with elements including:
#' * `trunc`: whether to truncate residual values below the argument "small"
#' * `small`: level of "small" values to guarentee numerical stability. If not specified, set dynamically based on the standard deviation of the outcome variable.
#' * `lambda`: For penalized regression, you can specify a level of lambda which will weight the penalty. If not set, will be determined based on 10-fold cross-validation.
#' * `output_quantiles`: whether to save fitted quantiles as part of the function output
#' * `calc_avg_me`: whether to return average marginal effects as part of the fitted object
#' * `lambda`: the penalization factor to be passed to penalized regression algorithms
#' @param std_err_control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`standard_errors`]. Possible arguments include:
#' * `se_method`: Method to use for standard errors, either "weighted_bootstrap",
#' "subsample", "bootstrap" or "custom" along with a specified subsampling method and
#' subsample percent. If specifying "custom", must also specify `subsampling_percent` and
#' `draw_weights`. If you specify "subsample", subsampling percent defaults to 0.2, but can be
#' changed. See details for details.
#' * `num_bs`: Number of bootstrap iterations to use, defaults to 100.
#' * `subsample_percent`: A number between 0 and one, specifying the percent of the data to subsample for standard error calculations
#' * `draw_weights`: Whether to use random exponential weights for bootstrap, either TRUE or FALSE
#' * `sampling_method` One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and bootstrapRows
#' samples with replacement.`
#' @param ... additional arguments, ignored for now
#' @details The qs function is a higher-level interface to fitting quantile spacings model,
#' handling both the quantile spacings regression, allowing the user to specify a number
#' of possible algorithms and methods for standard errors. It also supports
#' @importFrom assertthat assert_that
#' @importFrom SparseM model.matrix
#' @importFrom stats model.frame
#' @importFrom future nbrOfWorkers
#' @importFrom SparseM model.response
#' @export
qs <- function(formula, data = NULL,
               quantiles = c(0.9, 0.75, 0.5, 0.25, 0.1),
               baseline_quantile = 0.5,
               cluster_formula = NULL,
               weights = NULL,
               algorithm = "sfn",
               control = qs_control(),
               std_err_control = se_control(),
               parallel = TRUE,
               calc_se = TRUE,
               seed = NULL,
               ...) {


  algorithm <- check_algorithm(algorithm)

  assertthat::assert_that(length(baseline_quantile) == 1)

  m <- SparseM::model.matrix(formula, data)
  y <- SparseM::model.response(stats::model.frame(formula, data,
                                                  na.action = "na.pass"),
                               type = "numeric")

  if(is.null(control$small)) {
    control$small = pmax(sd(y)/5000, .Machine$double.eps)
  }

  quantiles <- unique(quantiles)

  if(!baseline_quantile %in% quantiles){
    quantiles <- c(quantiles, baseline_quantile)
  }

  reg_spec <- denseMatrixToSparse(m)
  reg_spec_var_names <- colnames(m)

  alpha <- sort(quantiles)

  jstar <- which(alpha == baseline_quantile)

  if(is.null(cluster_formula)) {
    cluster_matrix = NULL
  } else {
    cluster_matrix = SparseM::model.matrix(cluster_formula, data)
  }

  quantreg_fit <- quantreg_spacing(
    y = y,
    X = reg_spec,
    var_names = reg_spec_var_names,
    alpha = alpha,
    jstar = jstar,
    weights = weights,
    control = control,
    algorithm = algorithm,
    ...
  )

  if(grepl("lasso", algorithm)) {
    lambda = unlist(quantreg_fit$lambda)
  } else {
    lambda = NULL
  }

  if(calc_se) {
    if(!is.null(control$subsample_percent)) {
      assertthat::assert_that(control$subsample_percent > 0)
      assertthat::assert_that(control$subsample_percent <= 1)
    }
    se = standard_errors(
      y = y,
      X = reg_spec,
      cluster_matrix = cluster_matrix,
      var_names = reg_spec_var_names,
      weights = weights,
      alpha = alpha,
      jstar = jstar,
      algorithm = algorithm,
      std_err_control = std_err_control,
      parallel = parallel,
      control = qs_control(control$trunc, control$small, lambda = lambda,
                           output_quantiles = F, calc_avg_me = F),
      ...)
  } else {
    se = list(
      quant_cov = matrix(NA, nrow = ncol(reg_spec) * length(quantiles),
                         ncol = ncol(reg_spec) * length(quantiles))
    )
  }

  rv = list('quantreg_fit' = quantreg_fit,
            'se' = se,
            'specs' = list('formula' = formula,
                           'X' = data,
                           'Y' = y,
                           'alpha' = alpha,
                           'jstar' = jstar,
                           'trunc' = trunc,
                           'small' = control$small,
                           'subsample_percent' = std_err_control$subsample_percent,
                           'draw_weights' = std_err_control$draw_weights,
                           'num_bs' = std_err_control$num_bs,
                           'parallel' = parallel,
                           'coef_names' = colnames(m),
                           'algorithm' = algorithm))

  structure(rv, class = "qs")
}

#' creates a table of summary output for a qs object
#' @param object an object of class qs
#' @param ... other arguments to pass to summary
#' @importFrom utils head
#' @export
summary.qs = function(object, ...){

  quant_betas <- object$quantreg_fit$coef
  pseudo_r <- object$quantreg_fit$pseudo_r
  quant_se <- diag(object$se$quant_cov)^(0.5)

  coef_names <- object$specs$coef_names

  quant_out_betas <- t(matrix(quant_betas, ncol = length(object$specs$alpha)))
  quant_out_ses <- t(matrix(quant_se, ncol = length(object$specs$alpha)))

  baseline_quantile <- object$specs$jstar

  baseline_beta <- quant_out_betas[baseline_quantile,]
  baseline_se <- quant_out_ses[baseline_quantile,]

  baseline_mat <- data.frame(Variable = coef_names,
                             Quantile = object$specs$alpha[baseline_quantile],
                             Coefficient = unlist(baseline_beta), SE = baseline_se)

  quant_out_betas[baseline_quantile,] <- rep(0, ncol(quant_out_betas))
  quant_out_ses[baseline_quantile,] <- rep(NA, ncol(quant_out_ses))

  q_vec <- c()
  name_vec <- c()
  beta_vec <- c()
  se_vec <- c()

  for(i in 1:nrow(quant_out_betas)) {

    for(j in 1:ncol(quant_out_betas)) {

      id = i * ncol(quant_out_betas) - (ncol(quant_out_betas) - j)

      beta_vec[id] <- quant_out_betas[i, j][[1]]

      name_vec[id] <- coef_names[j]

      se_vec[id] <- quant_out_ses[i, j][[1]]

      q_vec[id] <- object$specs$alpha[i]
    }

  }

  quant_out_mat <- data.frame(Variable = name_vec, Quantile = q_vec,
                              Coefficient = beta_vec, `Standard Error` = se_vec)


  if(grepl("lasso", object$specs$algorithm)) {
    quant_out_mat <- quant_out_mat[which(quant_out_mat$Coefficient != 0),]
    baseline_mat <- baseline_mat[which(baseline_mat$Coefficient != 0),]
  }

  final_output <- list(
    baseline_coefs = baseline_mat,
    spacing_coefs = quant_out_mat[which(quant_out_mat$Quantile != object$specs$alpha[baseline_quantile]),],
    R2 = list(psuedo_r = pseudo_r),
    algorithm = object$specs$algorithm,
    quantiles = object$specs$alpha
  )

  structure(final_output, class = "qs_summary")
}

#' Capture print output
#' @param x object to capture
#' @param ... other argument to the print function
#' @importFrom testthat capture_output
capture_output <- function(x, ...) {
  testthat::capture_output(print(x, ...))
}

#' Round if x is numeric, otherwise don't
#' @param df data.frame whose columns I want to round
#' @param d number of digits
round_if <- function(df, d){
  rdf <- as.data.frame(lapply(df, FUN = function(x) {
    if(is.numeric(x)){
      signif(x, d)
    } else {
      x
    }
  }))
  structure(rdf, class = class(df))
}

#' Print qs summary
#' @param x fit from qs
#' @param digits number of digits to print
#' @param ... additional arguments, ignored for now
#' @export
print.qs <- function(x, digits = 4, ...) {

  d <- digits

  # really poor solution, going with it for now
  x <- summary(x, ...)

  spacing_coefs <- round_if(
    x$spacing_coefs, d
    )


  s_coefs <- capture_output(spacing_coefs)
  b_coefs <- capture_output(round_if(x$baseline_coefs, d))

  cat("Baseline Coefficients:\n",
      b_coefs, "\n\n")

  cat("Spacings Coefficients:\n",
      s_coefs, "\n\n")

  if(grepl("lasso", x$algorithm)) {
    cat("Only displaying non-zero coefficients.")
  }
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


#' Pad vector of strings based on longest length
#' @param x string vector
pad_strings <- function(x) {
  max_len <- max(nchar(x))

  for(i in 1:length(x)) {
    x[i] <- paste0(paste0(rep(" ", max_len - nchar(x[i])), collapse = ""), x[i])
  }

  x
}

#' Print qs summary
#' @param x fit from qs
#' @param digits number of digits to print
#' @param ... additional arguments, ignored for now
#' @export
print.qs_summary <- function(x, digits = 4, ...) {

  d <- digits

  spacing_coefs <- round_if(
    x$spacing_coefs, d
  )


  s_coefs <- capture_output(spacing_coefs)
  b_coefs <- capture_output(round_if(x$baseline_coefs, d))


  psuedo_r2 <- signif(x$R2$psuedo_r, d)
  psuedo_message <- paste0(paste0("    ",
                                  unique(x$quantiles),
                                      ":\t",
                                      psuedo_r2), collapse = "\n")

  cat("Baseline Quantile Coefficients:\n",
      b_coefs, "\n\n")

  cat("Quantile Spacing Coefficients:\n",
      s_coefs, "\n\n")

  cat(paste0("Quantile Psuedo-R2:\n",
       psuedo_message, "\n\n"))

}

#' Predict quantiles given fitted spacings model
#' @param object fitted `qs` model
#' @param newdata optional new data to pass
#' @param ... additional arguments, ignored for now
#' @return This function returns a set of predictive quantiles
#' for each data point, at the quantiles where spacings were
#' estimated.
#' @importFrom SparseM model.matrix
#' @importFrom stats as.formula
#' @importFrom stats delete.response
#' @importFrom stats terms
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats resid
#' @importFrom stats sd
#' @export
predict.qs <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) {
    if(!is.null(object$quantreg_fit$quantiles)) {
      p_q <- object$quantreg_fit$quantiles
    } else {
      X <- SparseM::model.matrix(stats::as.formula(object$specs$formula),
                                 data = object$specs$X)
      p_q <- spacings_to_quantiles(matrix(as.numeric(object$quantreg_fit$coef),
                                        ncol = length(object$specs$alpha)),
                                 X,
                                 jstar = object$specs$jstar)
    }
  } else {
    ff <- stats::as.formula(object$specs$formula)
    tt <- stats::terms(ff)
    tt <- stats::delete.response(tt)
    X <- stats::model.matrix(tt,
                              data = newdata)
    p_q <- spacings_to_quantiles(matrix(as.numeric(object$quantreg_fit$coef),
                                      ncol = length(object$specs$alpha)),
                               X,
                               jstar = object$specs$jstar)
  }

  colnames(p_q) <- object$specs$alpha
  p_q
}

#' Method for getting coefficients from fitted qs model
#' @param object fitted qs model
#' @param ... currently ignored
#' @export
coef.qs <- function(object, ...) {
  m <- t(matrix(unlist(object$quantreg_fit$coef), ncol = length(object$specs$alpha)))
  rownames(m) <- object$specs$alpha
  colnames(m) <- object$specs$coef_names
  m
}
