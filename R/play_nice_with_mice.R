
#' Imputation function to be used with the mice packe
#' @param y vector to be imputed
#' @param ry indicator for complete cases
#' @param x independent variables
#' @param wy cases to be imputed
#' @param ... other arguments to be passed to quantreg_spacing
mice.impute.qs <- function(y, ry, x, wy = NULL,
                           quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           baseline_quantile = 0.5,
                           algorithm = "rq.fit.lasso",
                           tails = "gaussian",
                           parallel = F,
                           calc_se = T,
                           weights = NULL,
                           control = qs_control(),
                           std_err_control = se_control(se_method = "weighted_bootstrap"),
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
      lambda = unlist(quantreg_fit$lambda)
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
                           output_quantiles = F, calc_avg_me = F),
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

  rng <- distributional_effects(out, tails = tails, alphas = alpha, parallel = parallel)
  rng$r(1)
}
