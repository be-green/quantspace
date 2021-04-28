#' Bin a variable along a range
#' @param x variable to create binned sequence along
#' @param size size of bin
#' @param trim percent of obs to trim from both ends of variable
#' @importFrom assertthat assert_that
#' @importFrom stats sd
#' @importFrom stats quantile
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

#' Computes means for various slices of a regression_spec by betas product.
#' @param x Regression specification.
#' @param betas Fitted quantile coefficients for specification.
#' @param alpha Quantiles targeted.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param trim Fraction (0 to 0.5) of observations to be trimmed from each end of
#' a given column in the regression_spec by betas product before the mean is computed.
#' @return Matrix with rows containing means for each quantile for various slices
#'  of the specification.
#' @export
avg_spacing = function(x, betas, alpha, jstar, trim=0){

  x_beta_prod <- as.matrix(x %*% betas)

  colnames(x_beta_prod) <- as.character(alpha)
  x_beta_prod[,-jstar] <- exp(x_beta_prod[,-jstar])

  apply(x_beta_prod, 2, mean, trim = trim)
}

#' Calculates the marginal effects  of an N x p matrix (wide-format) of qreg coefficients
#' @param qreg_coeffs wide-format of calculated spacings point estimates
#' @param avg_spacings average spacings matrix, which can be calculated from setting calculateAvgME = TRUE or using
#'                 the avg_spacing function directly
#' @param j_star the first quantile the user wishes to predict (usually the middle one)
#' @param calc_se boolean value, indicating whether the user wishes to calculate the marginal effect standard errors
#' @param qreg_vcv_vec variance-covariance matrix from point estimates, only necessary if calculating standard errors
#' @return list of values: avgME: calculated average marginal effects
#' avgME_se (user-specified): standard errors on the calculated marginal effects
get_marginal_effects = function(qreg_coeffs,
                                avg_spacings,
                                j_star,
                                calc_se = TRUE,
                                qreg_vcv_vec = NULL){

  N = dim(qreg_coeffs)[1]
  p = dim(qreg_coeffs)[2]
  avgME = matrix(NA, N, p)

  if(!missing(qreg_vcv_vec)) {
    avgME_se = avg_spacings*0
  } else {
    avgME_se = c()
  }

  # matrix with transformations of data that give marginal effects;
  # ME = R_matrix * qreg_coeffs
  R_matrix = array(0, dim = c(p,p,N))

  # calculating marginal effects here
  R_matrix[,j_star,] = 1

  for(jj in 1:N) {

    for(kk in 1:(j_star-1)){
      R_matrix[1:(j_star-kk),j_star-kk,jj] = -avg_spacings[j_star-kk]
    }
    for(kk in 1:(p-j_star)) {
      R_matrix[(j_star+kk):p,(j_star+kk),jj] = avg_spacings[j_star+kk]
    }

    avgME[jj,] = R_matrix[,,jj] %*% qreg_coeffs[jj,]

    # calculating standard errors here.
    # Can't see an easy way to avoid a double loop here
    if(calc_se && !missing(qreg_vcv_vec)){
      for(kk in 1:p){
        avgME_se[kk,jj] = sqrt(R_matrix[kk,,jj] %*%
                                 matrix(qreg_vcv_vec[,jj],p,p) %*%
                                 t(matrix(R_matrix[kk,,jj,drop=FALSE], 1, p)))
      }
    }
  }

  if(calc_se) return(list('avgME' = avgME, 'avgME_se' = avgME_se))
  else return(list('avgME' = avgME))
}

marginal_effects <- function(object, ...) {
  UseMethod("marginal_effects")
}

marginal_effects.qs <- function(fit, data, ...) {
  X <- stats::model.matrix(stats::as.formula(fit$specs$formula),
                           data = data)


}



#' Get marginal effects at a set of levels for the covariates
#' @param fit A fitted model from the `qs` function
#' @param data optional data.frame that specifies level of data to calculate
#' marginal effects
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

  me <- get_marginal_effects(reg_coefs, spacings, jstar)$avgME

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

