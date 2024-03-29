#' @rdname mice-helpers
#' @param data data to be interpolated by mice
#' @param quantiles Quantiles to target
#' @param ... other controls to be passed to the [`qs`] function
#' @importFrom stats complete.cases
#' @export
make_penalized_blots <- function(data,
                                 quantiles = c(0.01, 0.05, 0.1, 0.25, 0.5,
                                               0.75, 0.9, 0.95, 0.99),
                                 ...) {

  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  l <- vector("list", ncol(data))

  for (i in seq_along(l)) {
    l[[i]] <- alist()
  }

  names(l) <- colnames(data)

  use_qs <- lapply(data, function(x) is.numeric(x) & length(unique(x)) > 2)

  for(i in 1:length(l)) {
    if(use_qs[[i]]) {
      formula <- as.formula(paste0(colnames(data)[i], "~ ."))

      fit <- qs(formula, quantiles = quantiles,
                data = data[stats::complete.cases(data),],
                calc_se = F, algorithm = "lasso", ...)

      l[[i]]$control <- qs_control(lambda = unlist(fit$quantreg_fit$lambda),
                                   calc_avg_me = F, calc_r2 = F)
    }
  }
  l
}


#' Imputation function to be used with the mice package
#' @rdname mice-helpers
#' @param y vector to be imputed
#' @param ry indicator for complete cases
#' @param x independent variables
#' @param wy cases to be imputed
#' @param quantiles vector of quantiles to be estimated
#' @param baseline_quantile baseline quantile to measure spacings from (defaults to 0.5)
#' @param calc_se boolean, whether or not to calculate standard errors. Defaults to FALSE.
#' @param weights optional vector of weights for weighted quantile regression
#' @param parallel whether to run bootstrap in parallel
#' @param algorithm What algorithm to use for fitting underlying regressions.
#' Either one of "sfn", "br", "lasso", "post_lasso", or a function name which estimates
#' quantiles. Defaults to sfn for now.
#' @param tails what distribution to use when fitting the tails, either "gaussian" or "exponential"
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
#' the lower-level function called by [`standard_errors`]. Constructed via the [`se_control`] function.
#' Possible arguments include:
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
#' @param ... other arguments to be passed to quantreg_spacing
#' @examples
#' \dontrun{
#' library(mice)
#' x <- rnorm(10000)
#' x[sample(1:length(x), 100)] <- NA
#' x <- matrix(x, ncol = 10)
#'
#' # get optimal lambdas from CV search based on complete data
#' bl <- make_penalized_blots(x)
#'
#' # pass those to the lasso and get imputations
#' imputations = mice::mice(x, m = 10,
#'                          defaultMethod = c("qs", "logreg", "polyreg", "polr"),
#'                          blots = bl, algorithm = "lasso")
#'}
#' @export
#' @importFrom MASS mvrnorm
mice.impute.qs <- function(y, ry, x, wy = NULL,
                           quantiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                           baseline_quantile = 0.5,
                           algorithm = "two_pass",
                           tails = "gaussian",
                           parallel = F,
                           calc_se = F,
                           weights = NULL,
                           control = qs_control(calc_r2 = F,
                                                calc_avg_me = F),
                           std_err_control = se_control(),
                           ...) {
  if (is.null(wy)) {
    wy <- !ry
  }

  if(!baseline_quantile %in% quantiles){
    quantiles <- c(quantiles, baseline_quantile)
  }

  if(is.null(control$small)) {
    control$small = pmax(sd(y)/5000, .Machine$double.eps)
  }

  algorithm <- check_algorithm(algorithm)

  alpha <- sort(quantiles)
  jstar <- which(alpha == baseline_quantile)

  x <- cbind(1, x)

  fit <- quantreg_spacing(
    y = y[which(ry)],
    X = x[which(ry),],
    var_names = paste0("V", 1:ncol(x)),
    alpha = alpha,
    jstar = jstar,
    weights = weights,
    control = control,
    algorithm = algorithm,
    ...
  )

  if(calc_se) {

    if(grepl("lasso", algorithm)) {
      lambda = unlist(fit$lambda)
    } else {
      lambda = NULL
    }

    se = standard_errors(
      y = y[which(ry)],
      X = x[which(ry),],
      cluster_matrix = matrix(1, nrow = sum(ry)),
      var_names = paste0("V", 1:ncol(x)),
      weights = weights,
      alpha = alpha,
      jstar = jstar,
      algorithm = algorithm,
      std_err_control = std_err_control,
      parallel = parallel,
      control = qs_control(control$trunc, control$small, lambda = lambda,
                           output_quantiles = F, calc_avg_me = F, calc_r2 = F),
      ...)


    coefMat = MASS::mvrnorm(n = 1,
                            mu = unlist(fit$coef),
                            Sigma = se$quant_cov)

    coefMat <- matrix(coefMat, ncol = length(alpha))

  } else {
    coefMat <- matrix(unlist(fit$coef), ncol = length(alpha))
  }


  out <- spacings_to_quantiles(spacingCoef = coefMat,
                               data = x[which(wy),],
                               jstar = jstar)

  suppressWarnings({
    rng <- distributional_effects(out, tails = tails, alphas = alpha, parallel = parallel)
  })

  rng$r(1)
}
