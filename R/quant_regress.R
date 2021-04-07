#' Sparse Regression Quantile Fitting with Weights
#' @details A wrapper around the rq.fit.sfn function from the quantreg package,
#' extended to allow for  a user-supplied starting value and weights
#' @importFrom quantreg  rq.fit.sfn
#' @importFrom quantreg sfn.control
#' @importFrom quantreg sfnMessage
#' @param X structure of the design matrix X stored in csr format
#' @param y outcome vector
#' @param tau desired quantile
#' @param rhs the right-hand-side of the dual problem; regular users shouldn't need to specify this,
#' but in special cases can be quite usefully altered to meet special needs.
#' See e.g. Section 6.8 of Koenker (2005).
#' @param sv startin value for optimization, useful when bootstrapping
#' @param control control parameters for fitting routines: see [quantreg::sfn.control()]
#' @param weight_vec Optional vector of weights for regression
#' @param lambda ignored
#' @param ... other parameters, ignored
#' @export
rq.fit.sfn_start_val <- function(X,y,tau=.5,
                                 rhs = (1-tau)*c(t(a) %*% rep(1,length(y))),
                                 control,
                                 sv,
                                 weight_vec = NULL,
                                 lambda,
                                 ...) {
  a <- X
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
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weight_vec,`*`)
    a <- as(as.vector(weight_vec), "matrix.diag.csr") %*% a
  }

  u <- rep(1,length=n)
  x <- rep((1-tau),length=n)
  nnzdmax <- nnza <- a@ia[n+1]-1
  iwmax <- 7*m+3
  ao <- t(a)
  e <- ao %*% a
  nnzemax <- e@ia[m+1]-1
  ctrl <- quantreg::sfn.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  nsubmax <- ctrl$nsubmax
  tmpmax <- ctrl$tmpmax
  nnzlmax <- ctrl$nnzlmax
  if (is.null(ctrl$nsubmax)) nsubmax <- nnzemax
  if (is.null(ctrl$tmpmax)) tmpmax <- 6 * m
  if (is.null(ctrl$nnzlmax)) nnzlmax <- 4 * nnzdmax
  wwm <- vector("numeric",3*m)
  s <- u - x
  if(missing(sv)){
    b1 <- solve(e, ao %*% y, tmpmax=tmpmax,nnzlmax=nnzlmax,nsubmax=nsubmax)
  }
  else {
    # note: LDWS flipped the sign here, since formula above yields OLS coeff * -1
    b1 = -sv
  }

  r <- y - a %*% b1
  z <- ifelse(abs(r)<ctrl$small,(r*(r>0)+ctrl$small),r*(r>0))
  w <- z - r
  wwn <- matrix(0,n,14)
  wwn[,1] <- r
  wwn[,2] <- z
  wwn[,3] <- w
  fit <- .Fortran("srqfn",
                  n = as.integer(n),
                  m = as.integer(m),
                  nnza = as.integer(nnza),
                  a = as.double(a@ra),
                  ja = as.integer(a@ja),
                  ia = as.integer(a@ia),
                  ao = as.double(ao@ra),
                  jao = as.integer(ao@ja),
                  iao = as.integer(ao@ia),
                  nnzdmax = as.integer(nnzdmax),
                  d = double(nnzdmax),
                  jd = integer(nnzdmax),
                  id = integer(m+1),
                  dsub = double(nnzemax+1),
                  jdsub = integer(nnzemax+1),
                  nnzemax = as.integer(nnzemax),
                  e = as.double(e@ra),
                  je = as.integer(e@ja),
                  ie = as.integer(e@ia),
                  nsumax = as.integer(nsubmax),
                  lindx = integer(nsubmax),
                  xlindx = integer(m+1),
                  nnzlmax = as.integer(nnzlmax),
                  lnz = double(nnzlmax),
                  xlnz = integer(m+1),
                  iw = integer(m*5),
                  iwmax = as.integer(iwmax),
                  iwork = integer(iwmax),
                  xsuper = integer(m+1),
                  tmpmax = as.integer(tmpmax),
                  tmpvec = double(tmpmax),
                  wwm = as.double(wwm),
                  wwn = as.double(wwn),
                  cachsz = as.integer(ctrl$cachsz),
                  level = as.integer( 8 ),
                  x = as.double(x),
                  s = as.double(s),
                  u = as.double(u),
                  c = as.double(y),
                  sol = as.double(b1),
                  rhs = as.double(rhs),
                  small = as.double(ctrl$small),
                  ierr = integer(1),
                  maxiter = as.integer(ctrl$maxiter),
                  time = double(7),
                  PACKAGE = "quantreg")[c("sol","ierr",
                                          "maxiter","time")]
  ierr <- fit$ierr
  if(!(ierr==0) && ctrl$warn.mesg)
    warning(quantreg::sfnMessage(ierr))
  coefficients <- -fit$sol

  if(any(is.na(coefficients) | is.nan(coefficients))) {
    stop("Solution did not converge, coefficients are NA.")
  }

  residuals <- -y - a %*% coefficients
  if (!is.null(weight_vec)){
    residuals <- residuals / weight_vec
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = ctrl,
       ierr = ierr,
       it = fit$maxiter,
       weight_vec = weight_vec,
       out = NA)
}

#' Version that complies with more general requirements
#' @param X structure of the design matrix X stored in csr format
#' @param y response vector
#' @param tau target quantile
#' @param weight_vec optional vector of weights
#' @param rhs right hand size of dual problem
#' @param control control parameters for fitting routines
#' @param lambda ignored
#' @param ... other arguments, ignored
#' @importFrom quantreg rq.fit.sfn
rq.fit.sfn <- function(X, y, tau = 0.5,
                       weight_vec = NULL,
                       rhs = (1-tau)*c(t(X) %*% rep(1,length(y))),
                       control,
                       lambda,
                       ...) {

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){

    n = nrow(X)
    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weight_vec
    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weight_vec,`*`)
    X <- as(as.vector(weight_vec), "matrix.diag.csr") %*% X

  }
  c(quantreg::rq.fit.sfn(a = X, y, tau, rhs, control), out = NA)
}

#' Version that complies with more general requirements
#' @param X structure of the design matrix X stored in csr format
#' @param y response vector
#' @param tau target quantile
#' @param weight_vec optional vector of weights
#' @param control ignored
#' @param lambda ignored
#' @param ... additional quantities passed to rq.fit.br
#' @importFrom quantreg rq.fit.br
#' @importFrom stats coef
#' @importFrom stats resid
rq.fit.br <- function(X, y, tau = 0.5,
                      weight_vec = NULL, control,
                      lambda, ...) {

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){

    n = nrow(X)
    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weight_vec,`*`)
    X <- as(as.vector(weight_vec), "matrix.diag.csr") %*% X
  }

  fit <- quantreg::rq.fit.br(x = X, y, tau, ...)

  list(coefficients = coef(fit),
       residuals = resid(fit),
       control = list(),
       ierr = 0,
       it = 0,
       weight_vec = weight_vec,
       out = NA)
}

#' Quantile Regression w/ Lasso Penalty
#' @param X Design matrix, X
#' @param y outcome variable, y
#' @param tau quantile to estimate
#' @param lambda penalty parameter
#' @param weight_vec optional vector of weights
#' @param scale_x whether to scale the design matrix before estimation
#' @param method method argument to be passed to [quantreg::rq]
#' @param ... other arguments to pass to [rqPen::rq.lasso.fit]
#' @importFrom rqPen rq.lasso.fit
#' @importFrom rqPen cv.rq.pen
rq.fit.lasso <- function(X, y, tau, lambda, weight_vec,
                         scale_x = T, method = "br", ...) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if(nrow(X) < 100) {
    nfold = round(nrow(X)/10)
  } else {
    nfold = 10
  }


  intercept <- get_intercept(X)

  if(ncol(X) - (intercept > 0) <= 1) {
    stop("Model must have more than 1 covariate when using lasso.")
  }


  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){

    n = nrow(X)

    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    X <- diag(weight_vec) %*% X
  }

  # get the intercept, and then scale everything except the intercept
  if (scale_x) {
    unscaled_X <- X
    X <- scale_for_lasso(X, intercept)
    mu_x <- attr(X, "scaled:center")
    sigma_x <- attr(X, "scaled:scale")
    attr(X, "scaled:center") <- NULL
    attr(X, "scaled:scale") <- NULL
  }

  if(is.null(lambda)) {
    message("No lambda provided, selecting based on 10-fold cross-validation.")
    suppressWarnings({
      cv_fit <- rqPen::cv.rq.pen(x = X[,-intercept],
                              y = y, tau = tau,
                              intercept = T,
                              penalty = "LASSO",
                              criteria = "CV",
                              nfolds = nfold)
    })
    lambda = cv_fit$lambda.min
  }

  est <- rqPen::rq.lasso.fit(x = X[,-intercept],
                             y = y,
                             tau = tau,
                             lambda = lambda,
                             intercept = T,
                             scalex = F,
                             method = method,
                             ...)

  coefficients <- est$coefficients
  coefficients <- reorder_coefficients(coefficients, intercept)

  if(scale_x) {
    coefficients <- rescale_coefficients(coefficients, mu_x, sigma_x, intercept)
    X = unscaled_X
  }

  residuals = y - X %*% coefficients

  if (!is.null(weight_vec)){
    residuals <- residuals / weight_vec
  }

  list(coefficients = coef(est),
       residuals = residuals,
       control = NA,
       ierr = 0,
       it = est$it,
       weight_vec = weight_vec,
       out = list(lambda = lambda))
}

#' Quantile Regression w/ Lasso Penalty
#' @param X Design matrix, X
#' @param y outcome variable, y
#' @param tau quantile to estimate
#' @param lambda penalty parameter
#' @param weight_vec optional vector of weights
#' @param scale_x whether to scale the design matrix before estimation
#' @param method method argument to be passed to [quantreg::rq]
#' @param ... other arguments to pass to [rqPen::rq.lasso.fit]
#' @importFrom rqPen rq.lasso.fit
rq.fit.post_lasso <- function(X, y, tau, lambda, weight_vec,
                         scale_x = T, method = "br", ...) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }

  intercept <- get_intercept(X)

  if(nrow(X) < 100) {
    nfold = round(nrow(X)/10)
  } else {
    nfold = 10
  }


  if(ncol(X) - (intercept > 0) <= 1) {
    stop("Model must have more than 1 covariate when using lasso.")
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){
    message("No lambda provided, selecting based on 10-fold cross-validation.")
    n = nrow(X)

    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    X <- diag(weight_vec) %*% X
  }

  # get the intercept, and then scale everything except the intercept
  if (scale_x) {
    unscaled_X <- X
    X <- scale_for_lasso(X, intercept)
    mu_x <- attr(X, "scaled:center")
    sigma_x <- attr(X, "scaled:scale")
    attr(X, "scaled:center") <- NULL
    attr(X, "scaled:scale") <- NULL
  }

  if(is.null(lambda)) {
    message("No lambda provided, selecting based on 10-fold cross-validation.")
    suppressWarnings({
      cv_fit <- rqPen::cv.rq.pen(x = X[,-intercept],
                                 y = y, tau = tau,
                                 intercept = T,nfolds = nfold,
                                 penalty = "LASSO",
                                 criteria = "CV")
    })
    lambda = cv_fit$lambda.min
  }

  est <- rqPen::rq.lasso.fit(x = X[,-intercept],
                             y = y,
                             tau = tau,
                             lambda = lambda,
                             intercept = T,
                             scalex = F,
                             method = method,
                             ...)

  coefficients <- est$coefficients
  coefficients <- reorder_coefficients(coefficients, intercept)

  if(scale_x) {
    coefficients <- rescale_coefficients(coefficients, mu_x, sigma_x, intercept)
    X = unscaled_X
  }

  not_zero = which(coefficients != 0)
  new_X = X[,not_zero]

  est <- quantreg::rq.fit(x = new_X, y = y, tau = tau, method = method)
  coefficients[not_zero] <- est$coefficients

  residuals = y - X %*% coefficients

  if (!is.null(weight_vec)){
    residuals <- residuals / weight_vec
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = NA,
       ierr = 0,
       it = est$it,
       weight_vec = weight_vec,
       out = list(lambda = lambda))
}




#' Estimate a single quantile regression
#' @param X specification matrix for X variables
#' @param y outcome variable
#' @param tau quantile to regress
#' @param algorithm which algorithm to use
#' @param ... other arguments to be passed to the algorithm
fitQuantileRegression <- function(X, y, tau, algorithm = "rq.fit.sfn_start_val", ...) {
  do.call(algorithm, args = list(
    X = X, y = y, tau = tau, ...
  ))
}

#' Runs quantile regression on residuals of the model (calculates spaces around jstar quantile)
#' @param reg_spec_data result of ensureSpecRank function; regression matrix with full rank
#' @param ehat current residuals; subset of which to be used as dependent column
#' @param ind_hat column vector indicating which rows to be used in quantile regression
#' @param tau estimated quantile
#' @param trunc Boolean value; if true, replace those dependent values less than small with small itself;
#' else, only use rows with residuals greater than small
#' @param small Value used with trunc; values less than small 'blow up' too greatly when logged
#' @param algorithm The name of a function which will estimate a quantile regression.
#' Defaults to rq.fit.sfn_start_val. Must be a string, as it is passed to `do.call`
#' @param weight_vec vector of optional weights
#' @param ... other arguments to the function specified by the algorithm argument
#' @import SparseM
#' @return List of estimated coefficients, warnings, iterations, and controls as in
#' standard quantile regression function
#' @export
regressResiduals = function(reg_spec_data,
                        ehat,
                        ind_hat,
                        tau,
                        trunc,
                        small,
                        weight_vec,
                        algorithm = "rq.fit.sfn_start_val",
                        ...) {

  resids = ehat[ind_hat]

  if (!is.null(weight_vec)){
    weight_vec = as.matrix(weight_vec)
  }

  if(trunc) {
    resids <- log(pmax(resids,small))
    spec_mat <- reg_spec_data$spec_matrix
  } else {
    weight_vec = weight_vec[resids > small]
    # if the dimensions for the regression match, we don't subset the
    # x matrix, otherwise we do
    if(nrow(reg_spec_data$spec_matrix) == sum(resids > small)) {
      spec_mat <- reg_spec_data$spec_matrix
    } else {
      spec_mat <- reg_spec_data$spec_matrix[resids > small,]
    }
    resids <- log(resids[resids > small])
  }

  j_model <- fitQuantileRegression(
    X = spec_mat,
    y = resids,
    tau = tau,
    algorithm = algorithm,
    weight_vec = as.vector(weight_vec),
    ...)

  return(j_model)
}


#' Computes coefficients for the quantile regression spacing method.
#' @param dep_col Column of response variable.
#' @param data Regression specification matrix.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param small Minimum size of residuals for computational accuracy.
#' @param trunc Boolean value; if true, replace those dependent values less than small with small itself;
#' else, only use rows with residuals greater than small
#' @param algorithm The name of a function which will estimate a quantile regression.
#' Defaults to rq.fit.sfn_start_val. Must be a string, as it is passed to `do.call`
#' @param start_list Starting values for regression optimization.
#' @param weight_vec vector of optional weights
#' @param outputQuantiles TRUE or FALSE, whether to output quantiles
#' @param calculateAvgME TRUE or FALSE, whether to output average marginal effects
#' @param lambda optional penalty parameter, ignored except for penalized regression
#' algorithms
#' @param ... other parameters passed to the algorithm
#' @import SparseM
#' @return
#' Returns a list of coefficients.
#' num_betas is an x by p matrix of estimated parameters for each supplied quantiles.
#' pseudo_r is a  1 by p matrix of psuedo R^2 values for each quantile estimate.
#' warnings is a 1 by p matrix of warnings produced by each quantile regression call.
#' iter: is a 1 by p matrix of iterations ran by each quantile regression call.
#' @export
quant_reg_spacing = function(
  dep_col,
  data,
  var_names,
  alpha,
  jstar,
  algorithm = "rq.fit.sfn_start_val",
  small = 1e-3,
  trunc = FALSE,
  start_list = NA,
  weight_vec = NULL,
  outputQuantiles = FALSE,
  calculateAvgME = FALSE,
  lambda = NULL,
  ...) {

  # number of observations that will default to an NA estimate of
  # a quantile.
  obs_thresh = 2

  if(nrow(data) < obs_thresh) {
    stop("Data contains fewer than ", obs_thresh, "observations.")
  }

  width = dim(data)[2]
  tau = alpha[jstar]
  p = length(alpha)

  #create logs for output
  count_log = list()
  length(count_log) = p
  warnings_log = list()
  length(warnings_log) = p
  iter_log = list()
  length(iter_log) = p
  pseudo_r = list()
  length(pseudo_r) = p
  model = list() # Collect the quantile regressions into a list
  length(model) = p

  out = as.list(rep(NA, length(alpha)))

  # check to see if regression matrix is sparse. If not, then turn into CSR matrix
  if(!is(data, 'matrix.csr')) {
    data = denseMatrixToSparse(data)
  }

  if(missing(var_names) | is.null(var_names)) {
    var_names <- paste0("v", 1:ncol(data))
  }

  tmpmax <- floor(1e5 + exp(-12.1)*(pmin(data@ia[width+1], max(data@ia), na.rm = T)-1)^2.35)

  # Ensure matrix is not rank deficient
  reg_spec_starting_data <- ensureSpecFullRank(spec_mat = data, col_names = var_names)



  if(!is.null(lambda)) {
    if(length(lambda) > 1) {
      star_lambda = lambda[jstar]
    } else {
      star_lambda = lambda
    }
  } else {
    star_lambda = NULL
  }

  # Calculate initial fit
  if(any(is.na(start_list))){ # if user supplied starting values, then use them
    star_model = fitQuantileRegression(
      X = reg_spec_starting_data$spec_matrix,
      y = dep_col,
      tau = tau,
      control = list(tmpmax = tmpmax),
      weight_vec = weight_vec,
      algorithm = algorithm,
      lambda = star_lambda,
      ...)
   } else {
     col_nums = getColNums(start_list, reg_spec_starting_data, alpha, jstar)
     sv = as.numeric(start_list[col_nums])
     star_model = fitQuantileRegression(
       X = reg_spec_starting_data$spec_matrix,
       y = dep_col,
       tau = tau,
       control = list(tmpmax = tmpmax),
       sv = as.numeric(start_list[col_nums]),
       weight_vec = weight_vec,
       algorithm = algorithm,
       lambda = star_lambda,
       ...
      )
   }

  printWarnings(star_model)

  ehat0 = star_model$residuals

  #Calculate R^2
  V <- sum(rho(u = ehat0, tau = tau, weight_vec = weight_vec))
  V0 <- fitQuantileRegression(X = as.matrix.csr(rep(1, length(dep_col))),
                             y = dep_col,
                             tau = tau,
                             weight_vec = weight_vec,
                             algorithm = 'rq.fit.sfn_start_val')$residuals
  V0 <- sum(rho(u = V0, tau = tau,weight_vec = weight_vec))

  #set column names
  coef_df <- as.data.frame(t(star_model$coefficients))
  colnames(coef_df) <- reg_spec_starting_data$var_names
  coef_df <- addMissingSpecColumns(
    coef_df,
    var_names)
  colnames(coef_df) <- paste(alpha[jstar], colnames(coef_df), sep="_")

  #log output for return
  pseudo_r[[jstar]] = (1 - V/V0)
  model[[jstar]] = coef_df
  warnings_log[[jstar]] = star_model$ierr
  iter_log[[jstar]] = star_model$it
  count_log[[jstar]] = dim(reg_spec_starting_data$spec_matrix)[1]
  out[[jstar]] = star_model$out

  # Estimate upper quantiles sequentially
  ehat = ehat0
  for (j in (jstar+1):p) {
    ind_hat = which(ehat > 0)

    if(!is.null(lambda)) {
      if(length(lambda) > 1) {
        j_lambda = lambda[j]
      } else {
        j_lambda = lambda
      }
    } else {
      j_lambda = NULL
    }

    if(length(ind_hat) <= obs_thresh) {
      warning("When estimating coefficients quantile ", alpha[j], " there were ", obs_thresh,
              " or fewer residuals found above previous quantile estimate. Returning NA.")

      # Update residuals
      coef = rep(NA, length(star_model$coefficients))
      coef_df <- as.data.frame(t(coef))

      #get column names
      colnames(coef_df) <- reg_spec_starting_data$var_names
      coef_df <- addMissingSpecColumns(
        coef_df,
        var_names)
      colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

      #log results
      model[[j]] = coef_df
      pseudo_r[[j]] = NA
      warnings_log[[j]] = NA
      iter_log[[j]] = NA
      count_log[[j]] = NA

      # Update residuals
      ehat = NA
    } else {

      # Determine quantile to estimate
      tau.t = (alpha[j] - alpha[j-1])/(1 - alpha[j-1])

      # Ensure the cut of the starting data that we take for
      # current spacing is not rank-deficient
      if(!trunc){
        reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[which(ehat > small),],
                                            col_names = reg_spec_starting_data$var_names)
      } else {
        # else, handle rank specification for typically-sized matrix
        reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[ind_hat,],
                                            col_names = reg_spec_starting_data$var_names)
      }

      # run quantile regression
      coef <- NULL
      if(!any(is.na(start_list))){ # if user specified a start model, then use it as input
        col_nums = getColNums(start_list, reg_spec_data, alpha, j)
        sv = as.numeric(start_list[col_nums])
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    sv = sv, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                                    small = small, control = list(tmpmax = tmpmax),
                                    weight_vec = weight_vec[ind_hat],
                                    algorithm = algorithm, lambda = j_lambda,
                                    ...)
      } else{
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    ind_hat = ind_hat, tau = tau.t, trunc =  trunc,
                                    small = small, control = list(tmpmax = tmpmax),
                                    weight_vec = weight_vec[ind_hat],
                                    algorithm = algorithm, lambda = j_lambda,
                                    ...)
      }
      printWarnings(j_model)

      #Calculate R^2
      V <- sum(rho(u = j_model$residuals, tau = tau.t, weight_vec = j_model$weight_vec))
      V0 <- regressResiduals(
        reg_spec_data = list(
          spec_matrix = as.matrix.csr(rep(1, length(ind_hat)))
        ),
        ehat = ehat, ind_hat = ind_hat,
        tau = tau.t, trunc = trunc,
        small = small, weight_vec = weight_vec[ind_hat],
        algorithm = 'rq.fit.sfn')
      V0 <- sum(rho(u = V0$residuals, tau = tau.t,
                    weight_vec = V0$weight_vec))

      # Update residuals
      coef = j_model$coefficients
      coef_df <- as.data.frame(t(coef))

      #get column names
      colnames(coef_df) <- reg_spec_data$var_names
      coef_df <- addMissingSpecColumns(
        coef_df,
        var_names)
      colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

      #log results
      model[[j]] = coef_df
      pseudo_r[[j]] = (1 - V/V0)
      warnings_log[[j]] = j_model$ierr
      iter_log[[j]] = j_model$it
      count_log[[j]] = dim(reg_spec_data$spec_matrix)[1]
      out[[j]] = j_model$out

      # Update residuals
      ehat = ehat - exp(
        as.matrix(
          data %*%
            unname(t(as.matrix(model[[j]])))))

    }

  }

  # Estimate lower quantiles sequentially
  ehat = ehat0
  for (j in (jstar-1):1) {

    if(!is.null(lambda)) {
      if(length(lambda) > 1) {
        j_lambda = lambda[j]
      } else {
        j_lambda = lambda
      }
    } else {
      j_lambda = NULL
    }

    # Ensure the cut of the starting data that we take for
    # current spacing is not rank-deficient
    # if not truncating, then ensure exact rows of the regression matrix is included
    # else, handle rank specification for typically-sized matrix
    if(!trunc) {
      ind_hat = which(-ehat > small)
    } else {
      ind_hat = which(ehat < 0)
    }

    if(length(ind_hat) <= obs_thresh) {
      warning("When estimating coefficients quantile ", alpha[j], " there were ", obs_thresh,
              " or fewer residuals found above previous quantile estimate. Returning NA.")
      # Update residuals
      coef = rep(NA, length(star_model$coefficients))
      coef_df <- as.data.frame(t(coef))

      #get column names
      colnames(coef_df) <- reg_spec_starting_data$var_names
      coef_df <- addMissingSpecColumns(
        coef_df,
        var_names)
      colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

      #log results
      model[[j]] = coef_df
      pseudo_r[[j]] = NA
      warnings_log[[j]] = NA
      iter_log[[j]] = NA
      count_log[[j]] = NA

      # Update residuals
      ehat = NA
    } else {

      # Determine quantile to estimate
      tau.t = (alpha[j + 1] - alpha[j])/(alpha[j + 1])

      reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[ind_hat,],
                                          col_names = reg_spec_starting_data$var_names)


      #run quantile regression
      if(!any(is.na(start_list))) { # if the user supplied starting values, then use them in the function
        col_nums = getColNums(start_list, reg_spec_data, alpha, j)
        sv = as.numeric(start_list[col_nums])
        j_model <- regressResiduals(reg_spec_data = reg_spec_data,
                                    ehat = -ehat,
                                    sv = sv,
                                    ind_hat = ind_hat,
                                    tau = tau.t, trunc = trunc, small = small,
                                    control = list(tmpmax = tmpmax),
                                    algorithm = algorithm, lambda = j_lambda,
                                    weight_vec = weight_vec[ind_hat], ...)
      } else {
        j_model <- regressResiduals(reg_spec_data = reg_spec_data,
                                    ehat = -ehat,
                                    ind_hat = ind_hat,
                                    tau = tau.t, trunc = trunc, small = small,
                                    control = list(tmpmax = tmpmax),
                                    algorithm = algorithm,
                                    weight_vec = weight_vec[ind_hat],
                                    lambda = j_lambda,
                                    ...)
      }

      printWarnings(j_model)

      #Calculate pseudo-R^2
      V <- sum(rho(u = j_model$residuals, tau = tau.t, weight_vec = j_model$weight_vec))
      V0 <- regressResiduals(reg_spec_data = list('spec_matrix' = as.matrix.csr(rep(1, length(ind_hat)))),
                             ehat = -ehat, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                             small = small, weight_vec = weight_vec[ind_hat], algorithm = 'rq.fit.sfn')
      V0 <- sum(rho(u = V0$residuals, tau = tau.t, weight_vec = V0$weight_vec))

      # Update residuals
      coef = j_model$coefficients
      coef_df <- as.data.frame(t(coef))

      #get column names
      colnames(coef_df) <- reg_spec_data$var_names
      coef_df <- addMissingSpecColumns(
        coef_df,
        var_names)
      colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

      #log results
      model[[j]] = coef_df
      pseudo_r[[j]] = (1 - V/V0)
      warnings_log[[j]] = j_model$ierr
      iter_log[[j]] = j_model$it
      count_log[[j]] = dim(reg_spec_data$spec_matrix)[1]
      out[[j]] = j_model$out

        ehat = ehat + exp(as.matrix(
          data %*%
            unname(t(as.matrix(model[[j]])))))
    }
  }
  rv = list('coef' = do.call(cbind, model),
            'pseudo_r' = do.call(cbind, pseudo_r),
            'warnings' = do.call(cbind, warnings_log),
            'iter' = do.call(cbind, iter_log),
            'counts' = do.call(cbind, count_log),
            'out' = out)

  # calculate average marginal effects if user-specified
  if(calculateAvgME){
    # calculate the average spacing (column-wise average)
    averageSpacing = calcAvgSpacing(data, rv$coef, alpha,
                                    jstar, trim=.01)
    qreg_coef = matrix(as.numeric(rv$coef), ncol = p)
    me = getMarginalEffects(qreg_coeffs = qreg_coef,
                              avg_spacings = as.numeric(averageSpacing[2,]),
                              j_star = jstar,
                              calcSE = FALSE)
    rv$me = t(as.data.frame(as.vector(me$avgME)))
    colnames(rv$me) <- colnames(rv$coef)
    rownames(rv$me) <- NULL
  }

  # calculate quantiles induced by spacings if user-specified
  if(outputQuantiles) {
    rv$quantiles = spacings_to_quantiles(spacingCoef = matrix(as.numeric(rv$coef),
                                                            ncol = p),
                                       data, jstar)
  }
  return(rv)
}

#' Get the data inside the s4 slot of this sparse matrix class
#' @param data some thing that could be a sparse matrix
get_underlying <- function(data) {
  if(is.matrix.csr(data)) {
    data@ra
  } else {
    data
  }
}

#' Compute quantiles given parameter coefficients and data
#' @param spacingCoef J by p matrix; row is number of variables, p is number of quantiles
#' @param data independent variables
#' @param jstar index of median quantiles
#' @return N by p matrix of quantiles
#' @export
spacings_to_quantiles <- function(spacingCoef, data, jstar) {

  if(any(is.na(get_underlying(data))) | any(is.na(spacingCoef))) {
    spacingCoef <- as.matrix(spacingCoef)
    data <- as.matrix(data)
  }

  p = dim(spacingCoef)[2]
  quantiles = matrix(NA, nrow = dim(data)[1], ncol = p)

  starResids = data %*% spacingCoef[,jstar]
  quantiles[,jstar] = starResids

  resids = starResids
  for (j in (jstar+1):p) {
    spacing = data %*% spacingCoef[,j]
    resids = resids + exp(spacing)
    quantiles[,j] = resids
  }

  resids = starResids
  for (j in (jstar-1):1) {
    spacing = data %*% spacingCoef[,j]
    resids = resids - exp(spacing)
    quantiles[,j] = resids
  }

  return(quantiles)
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
calcAvgSpacing = function(x, betas, alpha, jstar, trim=.01){

  num_quantiles <- length(alpha)
  num_betas <- length(betas) / num_quantiles
  betas_wide <- matrix(as.numeric(betas), num_betas, num_quantiles)

  x_beta_prod <- as.matrix(x %*% betas_wide)

  colnames(x_beta_prod) <- as.character(alpha)
  x_beta_prod[,-jstar] <- exp(x_beta_prod[,-jstar])

  beta_means <- apply(x_beta_prod, 2, mean)
  beta_trimmed_means <- apply(x_beta_prod, 2, FUN = function(k) {
    mean(k, trim = trim)
    })

  means_out <- rbind(beta_means, beta_trimmed_means)
  rownames(means_out) <- c("FitQ_or_S_mean_full_spec",
                           "FitQ_or_S_mean_full_spec_trimmed")
  return(means_out)
}

#' Calculates the marginal effects  of an N x p matrix (wide-format) of qreg coefficients
#' @param qreg_coeffs wide-format of calculated spacings point estimates
#' @param avg_spacings average spacings matrix, which can be calculated from setting calculateAvgME = TRUE or using
#'                 the calcAvgSpacing function directly
#' @param j_star the first quantile the user wishes to predict (usually the middle one)
#' @param calcSE boolean value, indicating whether the user wishes to calculate the marginal effect standard errors
#' @param qreg_vcv_vec variance-covariance matrix from point estimates, only necessary if calculating standard errors
#
#' @return list of values: avgME: calculated average marginal effects
#' avgME_se (user-specified): standard errors on the calculated marginal effects
#' @export
getMarginalEffects = function(qreg_coeffs,
                              avg_spacings,
                              j_star,
                              calcSE = TRUE,
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
    for(kk in 1:(j_star-1)) R_matrix[1:(j_star-kk),j_star-kk,jj] = -avg_spacings[j_star-kk]

    for(kk in 1:(p-j_star)) R_matrix[(j_star+kk):p,(j_star+kk),jj] = avg_spacings[j_star+kk]

    avgME[jj,] = R_matrix[,,jj] %*% qreg_coeffs[jj,]

    # calculating standard errors here.
    # Can't see an easy way to avoid a double loop here
    if(calcSE && !missing(qreg_vcv_vec)){
      for(kk in 1:p){
        avgME_se[kk,jj] = sqrt(R_matrix[kk,,jj] %*%
                                 matrix(qreg_vcv_vec[,jj],p,p) %*%
                                 t(matrix(R_matrix[kk,,jj,drop=FALSE], 1, p)))
      }
    }
  }

  if(calcSE) return(list('avgME' = avgME, 'avgME_se' = avgME_se))
  else return(list('avgME' = avgME))
}
