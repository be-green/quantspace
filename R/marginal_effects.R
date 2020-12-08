
#' Bin a variable along a range
#' @param x variable to create binned sequence along
#' @param size size of bin
#' @param trim percent of obs to trim from both ends of variable
#' @importFrom assertthat assert_that
bin_along_range <- function(x, size = NA, trim = 0.01) {
  assertthat::assert_that(is.numeric(x))
  assertthat::assert_that(all(!is.na(x)))

  max_x <- quantile(x, 1 - trim)
  min_x <- quantile(x, trim)

  # slightly arbitrary, but w/e
  if(is.na(size)) {
    size = sd(x)/3
  }

  seq(min_x, max_x, by = size)
}

marginal_effects <- function(object, ...) {
  UseMethod("me")
}

marginal_effects.qs <- function(fit, data, ...) {
  X <- stats::model.matrix(stats::as.formula(fit$specs$formula),
                           data = data)


}



#' Get marginal effects at a set of levels for the covariates
#' @param fit A fitted model from the `qs` function
#' @param data optional data.frame that specifies level of data to calculate
#' marginal effects
#' @param variable Which coefficient
#' @export
#' @details A simple function which returns
#' marginal effects for a given level of the dataset.
me <- function(fit, data) {
  X <- stats::model.matrix(stats::as.formula(fit$specs$formula),
                             data = data)

  jstar <- fit$specs$jstar

  reg_coefs <- t(coef(fit))
  spacings <- as.matrix(X %*% reg_coefs)
  spacings[,-jstar] <- exp(spacings[,-jstar])

  me <- getMarginalEffects(reg_coefs, spacings, jstar)$avgME

  colnames(me) <- fit$specs$alpha
  rownames(me) <- fit$specs$coef_names

  me
}

#' Get marginal effects at a set of levels for the covariates
#' @param fit A fitted model from the `qs` function
#' @param variable variable to calculate marginal effects over
#' @param data optional data.frame that specifies level of data to calculate
#' marginal effects
#' @param size What bin size to use when varying the variable of interest
#' @param trim What to trim the variable of interest at, 0 < trim < 0.5
#' @details This function defaults to using
#' the average defaults to average for all coefficients
#' except variable chosen. If you want to vary a covariate but keep everything
#' else fixed at a certain level, you can specify the data argument.
#' The size argument defaults to a bin size of 1/3 of the standard deviation
#' of the variable, and the trim defaults to using the 95th percentile
#' instead of the max because there may be large outliers. You can over-ride
#' by setting trim to 0, which will use the min and max.
#' @importFrom SparseM model.matrix
#' @importFrom stats as.formula
#' @export
me_by_variable <- function(fit, variable, data = NA, size = NA, trim = 0.05) {
  if(is.na(data)) {
    X = SparseM::model.matrix(stats::as.formula(fit$specs$formula),
                              data = fit$specs$X)
    data <- data.frame(t(colMeans(X)))
  }

  vardata <- bin_along_range(fit$specs$X[[variable]])

  data <- as.data.frame(do.call("rbind", lapply(rep(list(data[,-which(colnames(data) == variable)]),
                    length(vardata)), unlist)))

  data[[variable]] <- vardata

  var_me <- matrix(NA, nrow = length(vardata), ncol = length(fit$specs$alpha))
  for(i in 1:length(data[[variable]])) {
    this_level_me <- me(fit, data = data[i,])

    var_me[i,] <- this_level_me[which(rownames(this_level_me) == variable),]
  }

  var_me <- as.data.frame(cbind(vardata, var_me))

  colnames(var_me) <- c(variable, fit$specs$alpha)

  var_me
}

#' Get all marginal effects of variables in the fit
#' @param fit model fitted by `qs()`
#' @param variable which variable to calculate marginal effects on
#' @param data optional data.frame that specifies level of data to calculate
#' marginal effects
#' @param size What bin size to use when varying the variable of interest
#' @param trim What to trim the variable of interest at, 0 < trim < 0.5
#' @details This function defaults to using
#' the average defaults to average for all coefficients
#' except variable chosen. If you want to vary a covariate but keep everything
#' else fixed at a certain level, you can specify the data argument.
#' The size argument defaults to a bin size of 1/3 of the standard deviation
#' of the variable, and the trim defaults to using the 95th percentile
#' instead of the max because there may be large outliers. You can over-ride
#' by setting trim to 0, which will use the min and max.
#' By default, marginal effects will calculate marginal effects for
#' all variables.
marginal_effects <- function(fit, variable = "all",
                             data = NA, size = NA,
                             trim = 0.05) {

  if(variable == "all") {

    d = SparseM::model.matrix(head(fit$specs$X, 100))

    variable <- colnames(d)
  }

  all_me <- lapply(variable,
                   function(v) {
                     data.frame(v, me_by_variable(fit, v, data, size, trim))
                   })

  do.call("rbind", all_me)

}

