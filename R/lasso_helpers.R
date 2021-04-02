
#' Rescale coefficients estimated on scaled data to match unscaled data
#' @param coefs scaled coefficients
#' @param mu_x vector of means for all (non-intercept) X variables
#' @param sigma_x vector of scales for all (non-intercept) X variables
#' @param intercept an integer for which column in design matrix is the intercept
rescale_coefficients <- function (coefs, mu_x, sigma_x, intercept) {
  new_coefs <- c()
  new_slope_coefs <- c()
  if (intercept != 0) {
    int_coef <- coefs[intercept]
    slope_coefs <- coefs[-intercept]
    for (j in 1:length(slope_coefs)) {
      new_slope_coefs[j] <- slope_coefs[j]/sigma_x[j]
      int_coef <- int_coef - slope_coefs[j] * mu_x[j]/sigma_x[j]
    }
    new_coefs[intercept] <- int_coef
    new_coefs[setdiff(1:length(coefs), intercept)] <- new_slope_coefs
  } else {
    for (j in 1:length(coefs)) {
      new_coefs[j] <- coefs[j]/sigma_x[j]
    }
  }
  new_coefs
}

#' Finds which column of X has the intercept term
#' @param X design matrix for regression
get_intercept <- function(X) {
  which_cols <- 1:ncol(X)
  i = 1
  while(length(which_cols) > 1 & i <= nrow(X)) {
    which_cols <- which(X[i, which_cols] == 1)
    i = i + 1
  }

  if (length(which_cols) > 1) {
    stop("Two intercept columns specified in design matrix. Please remove one.")
  } else if (length(which_cols) == 0) {
    0
  } else {
    if(all(X[, which_cols] == 1)) {
      which_cols
    } else {
      0
    }
  }
}

#' Scale matrix for lasso regression
#' @param X design matrix
#' @param intercept column number for intercept
scale_for_lasso <- function(X, intercept) {
  if(intercept == 0) {
    scale(X)
  } else {
    scaled_X = scale(X[,setdiff(1:ncol(X), intercept)])

    out_X = matrix(nrow = nrow(X), ncol = ncol(X))
    out_X[,intercept] = rep(1, nrow(X))
    out_X[,setdiff(1:ncol(X), intercept)] = scaled_X

    # reattach attributes to the matrix that has the intercept added back
    at = attributes(scaled_X)
    at = subset(at, sapply(names(at), function(x) grepl(pattern = "scaled", x)))

    for(i in 1:length(at)) {
      attr(out_X, names(at)[i]) <- at[i][[1]]
    }
    out_X
  }
}

#' Reorder coefficients in case intercept wasn't the first term
#' @param coefficients vector of coefficients
#' @param intercept column number for intercept term
reorder_coefficients <- function(coefficients, intercept) {
  if(intercept == 0) {
    coefficients
  } else {
    reordered_coefs <- rep(NA, length(coefficients))
    reordered_coefs[-intercept] <- coefficients[-1]
    reordered_coefs[intercept] <- coefficients[1]
    reordered_coefs
  }
}
