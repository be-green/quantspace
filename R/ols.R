
#' Computes OLS coefficients without losing "sparse"-ness of inputs
#' @param a Regression specification matrix
#' @param y Dependent variable
#' @param weight_vec Optional vector of weights
#' @return OLS coefficients
ols_sparse_fit <- function(a, y, weight_vec = NULL) {

  y <- -y
  n <- length(y)
  m <- a@dimension[2]
  if(n != a@dimension[1]) {
    stop("Dimensions of design matrix and the response vector not compatible")
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){
    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * sqrt(weight_vec)

    # multiplying the a matrix by the weights
    a <- as(as.vector(sqrt(weight_vec)), "matrix.diag.csr") %*% a

  }

  nnzdmax <- nnza <- a@ia[n+1]-1
  iwmax <- 7*m+3
  ao <- t(a)
  e <- ao %*% a
  nnzemax <- e@ia[m+1]-1
  nsubmax <- nnzemax
  nnzlmax <- 4 * nnzdmax
  tmpmax <- floor(1e5 + exp(-12.1)*(a@ia[m+1]-1)^2.35)

  return(-1 * solve(e, ao %*% y, tmpmax=tmpmax, nnzlmax=nnzlmax, nsubmax=nsubmax))
}



#' Get R^2 of OLS fit
#' @param a Regression specification matrix
#' @param y Outcome vector
#' @param betas Parameters from OLS fit
#' @param weight_vec Optional vector of weights
get_ols_r_squared <- function(a, y, betas, weight_vec = NULL) {
  if (is.null(weight_vec)) {
    e <- y - (a %*% betas)
    mu <- mean(y)
  } else {
    # transforming regressors to incorporate weights
    e <- (y - (a %*% betas)) * sqrt(weight_vec)
    y <- y * sqrt(weight_vec)
    mu <- mean(y )
  }
  return(1 - ((t(e) %*% e) / (t(y - mu) %*% (y - mu))))
}
