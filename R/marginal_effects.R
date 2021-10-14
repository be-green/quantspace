#' Bin a variable along a range
#' @param x variable to create binned sequence along
#' @param size size of bin
#' @param trim percent of obs to trim from both ends of variable
#' @importFrom assertthat assert_that
#' @importFrom stats sd
#' @importFrom stats quantile
bin_along_range <- function(x, size = NA, trim = 0.01) {
  if(length(setdiff(x, c(0,1))) == 0) {
    return(c(0, 1))
  } else {
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
#' @importFrom stringr str_extract
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

  nms <- colnames(qreg_vcv_vec)
  if(is.null(nms)) {
    nms <- rownames(qreg_coeffs)
  }
  qtls <- stringr::str_extract(nms, "^[0-9\\.]+")
  coef_nms <- stringr::str_replace(nms, paste0(qtls, "_"), "")

  avgME_se = matrix(NA, N, p)

  # matrix with transformations of data that give marginal effects;
  # ME = R_matrix * qreg_coeffs
  R_matrix = array(0, dim = c(p,p,N))

  # calculating marginal effects here
  R_matrix[,j_star,] = 1

  for(jj in 1:N) {
    nm <- rownames(qreg_coeffs)[jj]
    covmat_rows = which(coef_nms == nm)

    for(kk in 1:(j_star-1)){
      R_matrix[1:(j_star-kk),j_star-kk,jj] = -avg_spacings[j_star-kk]
    }
    for(kk in 1:(p-j_star)) {
      R_matrix[(j_star+kk):p,(j_star+kk),jj] = avg_spacings[j_star+kk]
    }

    if(calc_se == F | length(covmat_rows) == p) {
      avgME[jj,] = R_matrix[,,jj] %*% qreg_coeffs[jj,]

      # calculating standard errors here.
      # Can't see an easy way to avoid a double loop here
      if(calc_se && !missing(qreg_vcv_vec)){
        for(kk in 1:p){
          avgME_se[jj,kk] = sqrt(R_matrix[kk,,jj] %*%
                                   qreg_vcv_vec[covmat_rows,
                                                covmat_rows] %*%
                                   matrix(R_matrix[kk,,jj], ncol = 1))
        }
      } else {
        for(kk in 1:p){
          avgME_se[jj,kk] = NA
        }
      }
    } else {
      avgME[jj,] = NA

      # calculating standard errors here.
      # Can't see an easy way to avoid a double loop here
      if(calc_se && !missing(qreg_vcv_vec)){
        for(kk in 1:p){
          avgME_se[jj,kk] = NA
        }
      }
    }

  }

    return(list('avgME' = avgME,
                'avgME_se' = avgME_se))
}

#' Get marginal effects at a set of levels for the covariates
#' @param fit A fitted model from the `qs` function
#' @param X model matrix that specifies level of data to calculate
#' marginal effects
#' @export
#' @details A simple function which returns
#' marginal effects for a given level of the dataset.
me <- function(fit, X) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }
  jstar <- fit$specs$jstar
  alphas <- fit$specs$alpha
  reg_coefs <- t(coef(fit))
  spacings <- avg_spacing(X, reg_coefs, alphas, jstar)

  if(any(is.na(fit$se$quant_cov))) {
    calcSE = F
  } else {
    calcSE = T
  }
  me <- get_marginal_effects(reg_coefs, spacings, jstar,
                             calc_se = calcSE,
                             qreg_vcv_vec = fit$se$quant_cov)

  # point estimates
  colnames(me[[1]]) <- fit$specs$alpha
  rownames(me[[1]]) <- fit$specs$coef_names

  # standard errors
  colnames(me[[2]]) <- fit$specs$alpha
  rownames(me[[2]]) <- fit$specs$coef_names

  me
}

#' Get all marginal effects of variables in the fit
#' @param fit model fitted by `qs()`
#' @param type one of "ame" (average marginal effects) or
#' "mea" (marginal effects at the average)
#' @param variable which variable to calculate marginal effects on
#' @param data optional data.frame that specifies level of data to calculate
#' marginal effects
#' @param trim What to trim the variable of interest at, 0 < trim < 0.5
#' @details The trim defaults to using the 95th percentile
#' instead of the max because there may be large outliers. You can over-ride
#' by setting trim to 0, which will use the min and max.
#' By default, marginal effects will calculate marginal effects for
#' all variables.
#' @importFrom stats terms
#' @export
marginal_effects <- function(fit,
                             type = "mea",
                             variable = "all",
                             data = NA,
                             trim = 0.05) {

  if(anyNA(data) & length(data) == 1) {
    data = stats::model.matrix(fit$specs$formula, data = fit$specs$X)
  }
  if(length(variable) == 1 && variable == "all") {
    variable <- setdiff(colnames(data), "(Intercept)")
  }


  if(type == "mea") {
    data = matrix(colMeans(data), nrow = 1)
  }
  all_me <- me(fit, X = data)
  if(length(variable) > 1 || variable != "all") {
    all_me$avgME <- all_me$avgME[which(rownames(all_me$avgME) %in% variable),, drop = F]
    all_me$avgME_se <- all_me$avgME_se[which(rownames(all_me$avgME_se) %in% variable),,  drop = F]
  }

  attr(all_me, "jstar") <- fit$specs$jstar

  # ff = stats::as.formula(fit$specs$formula)
  # attr(all_me, "outcome") <- all.vars(ff)[attr(stats::terms(ff, data = data),
  #                                                     "response")]

  attr(all_me, "type") <- type

  structure(all_me,
            class = "qs_me")

}

#' Print method for quantspace marginal effects
#' @param x marginal effects to print
#' @param ... other arguments for s3 consistency, ignored for now
#' @export
print.qs_me <- function(x, ...) {
  type = attr(x, "type")
  out_type <- switch(type,
                     ame = "Average Marginal Effects:\n",
                     mea = "Marginal Effects at the Average:\n")
  cat(out_type)
  # cat(paste0(names(x)[i],": \n"))
  mat <- matrix(nrow = nrow(x$avgME) + nrow(x$avgME_se),
                ncol = ncol(x$avgME))
  rn = c()
  sub_x_rn = rownames(x$avgME)
  for(i in 1:nrow(x$avgME)) {
      rn <- c(rn, c(sub_x_rn[i], paste(sub_x_rn[i], "SE")))
      mat[2 * i - 1,] <- x$avgME[i,]
      mat[2 * i,] <- x$avgME_se[i,]

  }
  rownames(mat) <- rn
  colnames(mat) <- colnames(x$avgME)
  print(mat, ...)


}
