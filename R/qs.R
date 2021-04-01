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
#' @param subsamplePct percent to subsample for standard error calculations
#' @param cluster_indices index of clusters for clustered standard errors
#' @param stratum_indices index of strata for clustered standard errors
#' @param draw_weights whether to use random exponential weights for bootstrap
#' @param num_bs number of bootstrap draws
#' @param parallel whether to run bootstrap in parallel
#' @param num_cores number of cores to use (defaults to setting from `getOption(mc.cores)`)
#' @param save_rows whether to save rows, see [subsampleStandardErrors]
#' @param seed what seed to use for replicable RNG
#' @param ... additional arguments, ignored for now
#' @importFrom assertthat assert_that
#' @importFrom SparseM model.matrix
#' @importFrom stats model.frame
#' @importFrom SparseM model.response
#' @export
qs <- function(formula, data = NULL,
               quantiles = c(0.9, 0.75, 0.5, 0.25, 0.1),
               baseline_quantile = 0.5,
               se_method = "boot",
               weight_vec = NULL,
               subsamplePct = 0.2,
               algorithm = "sfn",
               cluster_indices = NULL,
               stratum_indices = NULL,
               draw_weights = TRUE,
               num_bs = 100,
               parallel = TRUE,
               num_cores = getCores(),
               trunc = T,
               small = NULL,
               seed = NULL,
               ...) {

  if(!exists(algorithm)) {
    if(algorithm == "sfn") {
      algorithm = "rq.fit.sfn"
    } else if(algorithm == "lasso") {
      algorithm = "rq.fit.lasso"
    } else {
      stop(paste0("Algorithm not implemented in quantspace, and not a function" ,
                  " available in the current namespace."))
    }
  }

  assertthat::assert_that(length(baseline_quantile) == 1)

  if(is.null(num_cores)) {
    num_cores <- getCores()
  }


  m <- SparseM::model.matrix(formula, data)
  y <- SparseM::model.response(stats::model.frame(formula, data), type = "numeric")

  depCol <- y


  if(is.null(small)) {
    small = pmax(sd(y)/5000, .Machine$double.eps)
  }

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
    algorithm = algorithm,
    outputQuantiles = T,
    calculateAvgME = F,
    ...
  )

  if(subsamplePct > 1) subsamplePct = subsamplePct * 0.01

  if(se_method == "boot") {
    se = subsampleStandardErrors(
      dep_col = depCol,
      data = reg_spec,
      var_names = reg_spec_var_names,
      alpha = alpha,
      jstar = jstar,
      algorithm = algorithm,
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
      ...)
  } else {
    stop("This method is not yet implemented or integrated with qs.")
  }


  rv = list('quantreg_fit' = quantreg_fit,
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
                           'coef_names' = colnames(m),
                           'algorithm' = algorithm))

  structure(rv, class = "qs")
}

#' creates a table of summary output for a qs object
#' @param x an object of class qs
#' @param ... other arguments to pass to summary
#' @importFrom utils head
#' @export
summary.qs = function(x, ...){

  quant_betas <- x$quantreg_fit$coef
  pseudo_r <- x$quantreg_fit$pseudo_r
  quant_se <- diag(x$se$quant_cov)^(0.5)

  coef_names <- x$specs$coef_names

  quant_out_betas <- t(matrix(quant_betas, ncol = length(x$specs$alpha)))
  quant_out_ses <- t(matrix(quant_se, ncol = length(x$specs$alpha)))

  baseline_quantile <- x$specs$jstar

  baseline_beta <- quant_out_betas[baseline_quantile,]
  baseline_se <- quant_out_ses[baseline_quantile,]

  baseline_mat <- data.frame(Variable = coef_names,
                             Quantile = x$specs$alpha[baseline_quantile],
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

      q_vec[id] <- x$specs$alpha[i]
    }

  }

  quant_out_mat <- data.frame(Variable = name_vec, Quantile = q_vec,
                              Coefficient = beta_vec, `Standard Error` = se_vec)


  if(x$specs$algorithm == "rq.fit.lasso") {
    quant_out_mat <- quant_out_mat[which(quant_out_mat$Coefficient != 0),]
    baseline_mat <- baseline_mat[which(baseline_mat$Coefficient != 0),]
  }

  final_output <- list(
    baseline_coefs = baseline_mat,
    spacing_coefs = quant_out_mat[which(quant_out_mat$Quantile != x$specs$alpha[baseline_quantile]),],
    R2 = list(psuedo_r = pseudo_r),
    algorithm = x$specs$algorithm
  )

  structure(final_output, class = "qs_summary")
}

#' Capture print output
#' @param x object to capture
#' @importFrom testthat capture_output
capture_output <- function(x, ...) {
  testthat::capture_output(print(x, ...))
}

#' Round if x is numeric, otherwise don't
#' @param df data.frame whose columns I want to round
#' @param d number of digits
#' @importFrom purrr map_df
round_if <- function(df, d){
  rdf <- purrr::map_df(df, .f = function(x) {
    if(is.numeric(x)){
      signif(x, d)
    } else {
      x
    }
  })
  structure(rdf, class = class(df))
}

#' Print qs summary
#' @param x fit from qs
#' @param d number of digits to print
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

  if(x$algorithm == "rq.fit.lasso") {
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
#' @param d number of digits to print
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
                                  unique(x$spacing_coefs$Quantile),
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
#' @export
predict.qs <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) {
    if(!is.null(object$quantreg_fit$quantiles)) {
      p_q <- object$quantreg_fit$quantiles
    } else {
      X <- SparseM::model.matrix(stats::as.formula(object$specs$formula),
                                 data = object$specs$X)
      p_q <- spacingsToQuantiles(matrix(as.numeric(object$quantreg_fit$coef),
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
    p_q <- spacingsToQuantiles(matrix(as.numeric(object$quantreg_fit$coef),
                                      ncol = length(object$specs$alpha)),
                               X,
                               jstar = object$specs$jstar)
  }

  colnames(p_q) <- object$specs$alpha
  p_q
}

#' Method for getting coefficients from fitted qs model
#' @param x fitted qs model
#' @param ... currently ignored
#' @export
coef.qs <- function(x, ...) {
  m <- t(matrix(unlist(x$quantreg_fit$coef), ncol = length(x$specs$alpha)))
  rownames(m) <- x$specs$alpha
  colnames(m) <- x$specs$coef_names
  m
}
