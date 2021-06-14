#' Sparse Regression Quantile Fitting with Weights
#' @details A wrapper around the rq.fit.sfn function from the quantreg package,
#' extended to allow for  a user-supplied starting value and weights
#' @importFrom quantreg  rq.fit.sfn
#' @importFrom quantreg sfn.control
#' @importFrom quantreg sfnMessage
#' @importFrom SparseM as.matrix.csr
#' @importFrom SparseM as.matrix
#' @param X structure of the design matrix X stored in csr format
#' @param y outcome vector
#' @param tau desired quantile
#' @param rhs the right-hand-side of the dual problem; regular users shouldn't need to specify this,
#' but in special cases can be quite usefully altered to meet special needs.
#' See e.g. Section 6.8 of Koenker (2005).
#' @param sv starting value for optimization, useful when bootstrapping
#' @param control control parameters for fitting routines: see [quantreg::sfn.control()]
#' @param weights Optional vector of weights for regression
#' @param lambda ignored
#' @param ... other parameters, ignored
#' @export
rq.fit.sfn_start_val <- function(X,y,tau=.5,
                                 rhs = (1-tau)*c(t(a) %*% rep(1,length(y))),
                                 control,
                                 sv,
                                 weights = NULL,
                                 lambda,
                                 ...) {
  if(inherits(X, "matrix")) {
    X <- denseMatrixToSparse(X)
  }

  a <- X
  y <- -y
  n <- length(y)
  m <- a@dimension[2]
  if(n != a@dimension[1]) {
    stop("Dimensions of design matrix and the response vector not compatible")
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){
    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * weights

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weights,`*`)
    a <- SparseM::as.matrix.csr(diag(as.vector(weights))) %*% a
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
    b1 <- solve(SparseM::as.matrix(e), ao %*% y, tmpmax=tmpmax,nnzlmax=nnzlmax,nsubmax=nsubmax)
  }
  else {
    # note: LDWS flipped the sign here, since formula above yields OLS coeff * -1
    b1 = -sv
  }
  r <- y - SparseM::as.matrix(a %*% b1)
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
  if (!is.null(weights)){
    residuals <- residuals / weights
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = ctrl,
       ierr = ierr,
       it = fit$maxiter,
       weights = weights)
}

#' Version that complies with more general requirements
#' @param X structure of the design matrix X stored in csr format
#' @param y response vector
#' @param tau target quantile
#' @param weights optional vector of weights
#' @param rhs right hand size of dual problem
#' @param control control parameters for fitting routines
#' @param lambda ignored
#' @param ... other arguments, ignored
#' @importFrom quantreg rq.fit.sfn
#' @importFrom SparseM t
#' @export
rq.fit.sfn <- function(X, y, tau = 0.5,
                       weights = NULL,
                       rhs = (1-tau)*c(t(X) %*% rep(1,length(y))),
                       control,
                       lambda,
                       ...) {

  if(inherits(X, "matrix")) {
    X <- denseMatrixToSparse(X)
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){

    n = nrow(X)
    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weights
    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weights,`*`)
    X <- SparseM::as.matrix.csr(diag(as.vector(weights))) %*% X

  }
  quantreg::rq.fit.sfn(a = X,y, tau, rhs, control)
}

#' Version that complies with more general requirements
#' @param X structure of the design matrix X stored in csr format
#' @param y response vector
#' @param tau target quantile
#' @param weights optional vector of weights
#' @param control ignored
#' @param lambda ignored
#' @param ... additional quantities passed to rq.fit.br
#' @importFrom quantreg rq.fit.br
#' @importFrom stats coef
#' @importFrom stats resid
#' @export
rq.fit.br <- function(X, y, tau = 0.5,
                      weights = NULL, control,
                      lambda, ...) {

  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){

    n = nrow(X)
    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * weights

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weights,`*`)
    X <- as(as.vector(weights), "matrix.diag.csr") %*% X
  }

  fit <- quantreg::rq.fit.br(x = X, y, tau, ...)

  list(coefficients = coef(fit),
       residuals = resid(fit),
       control = list(),
       ierr = 0,
       it = 0,
       weights = weights)
}

#' Quantile Regression approximated w/ huber loss
#' @param X design matrix
#' @param y outcome vector
#' @param tau target quantile
#' @param weights optional weight vector
#' @param control ignored for now
#' @param lambda ignored for now
#' @param smoothing_window neighborhood around 0 which is
#' smoothed by either typical least squares or appropriately
#' tilted least squares loss function
#' @param beta_tol stopping rule based on max value of gradient
#' @param check_tol stopping rule based on change in the loss function
#' @param maxiter largest number of iterations allowed
#' @param n_samples number of observations to use in "warmup" regression
#' @param init_beta initial guess at betas
#' @param scale whether to scale x and y variables in regression
#' @param intercept optional integer indicating intercept column
#' that identifies initial values
#' @param ... other arguments, ignored for now
#' @export
#' @importFrom stats rnorm
rq.fit.agd <- function(X, y, tau = 0.5,
                       weights = NULL, control,
                       lambda, smoothing_window = 1e-16,
                       beta_tol = 1e-16,
                       check_tol = 1e-16,
                       maxiter = 10000,
                       n_samples = min(c(ceiling(nrow(X)/10),
                                         10000)),
                       init_beta = NULL,
                       scale = 1,
                       intercept = NULL,
                       ...) {


  if(!inherits(X, "matrix")) {
    X <- as.matrix(X)
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){

    n = nrow(X)
    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * weights

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weights,`*`)
    X <- weights * X
  }
  if(is.null(intercept)) {
    intercept <- get_intercept(X)
  }
  if(is.null(init_beta)) {
    init_beta = stats::rnorm(ncol(X))
  }

  if(n_samples < 10) {
    n_samples = min(10, length(y))
  }

  samples = sample(1:length(y), n_samples)



  fit = fit_approx_quantile_model(X = X, y = y,
                          X_sub = X[samples,],
                          y_sub = y[samples],
                          tau = tau,
                          mu = smoothing_window, init_beta = init_beta,
                          maxiter = maxiter,
                          beta_tol = beta_tol,check_tol = check_tol,
                          intercept = intercept, num_samples = n_samples,
                          scale = scale,
                          )

  list(coefficients = fit,
       residuals = y - X %*% fit,
       control = list(),
       ierr = 0,
       it = 0,
       weights = weights)
}

#' Quantile Regression w/ Lasso Penalty
#' @param X Design matrix, X
#' @param y outcome variable, y
#' @param tau quantile to estimate
#' @param lambda penalty parameter
#' @param weights optional vector of weights
#' @param scale_x whether to scale the design matrix before estimation
#' @param method method to use when fitting underlying quantile regression algorithm
#' @param nfold number of folds to use when cross-validating
#' @param parallel whether to run cv search in parallel, if applicable
#' @param ... other arguments to pass to underlying fitting algorithm
#' @param nlambda number of lambdas to search over.
#' @export
rq.fit.lasso <- function(X, y, tau, lambda, weights,
                         scale_x = T, method = "agd", nfold = 10,
                         nlambda = 50, parallel = F, ...) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }

  intercept <- get_intercept(X)

  if(ncol(X) - (intercept > 0) <= 1) {
    stop("Model must have more than 1 covariate when using lasso.")
  }

  if(nrow(X) < nfold) {
    nfold = nrow(X)
  }


  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){

    n = nrow(X)

    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weights

    # pre-multiplying the a matrix by a diagonal matrix of weights
    X <- diag(weights) %*% X
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
    message("No lambda provided, selecting based on ",nfold,"-fold cross-validation.")
    suppressWarnings({
      cv_fit <- lasso_cv_search(x = X[,-intercept], y = y , tau = tau,
                                method = method,
                                intercept = TRUE, nfolds = nfold,
                                foldid = NULL, nlambda = nlambda,
                                eps = 1e-04, init.lambda = 1,
                                parallel = parallel,
                                scalex = FALSE,
                                ...)
    })
    lambda = cv_fit$lambda.min
  }

  est <- fit_lasso(x = X[,-intercept],
                             y = y,
                             tau = tau,
                             lambda = lambda,
                             intercept = TRUE,
                             scalex = FALSE,
                             method = method,
                             ...)

  coefficients <- est$coefficients
  coefficients <- reorder_coefficients(coefficients, intercept)

  if(scale_x) {
    coefficients <- rescale_coefficients(coefficients, mu_x, sigma_x, intercept)
    X = unscaled_X
  }

  residuals = y - X %*% coefficients

  if (!is.null(weights)){
    residuals <- residuals / weights
  }

  list(coefficients = coef(est),
       residuals = residuals,
       control = NA,
       ierr = 0,
       it = est$it,
       weights = weights,
       lambda = lambda)
}

#' Quantile Regression w/ Lasso Penalty
#' @param X Design matrix, X
#' @param y outcome variable, y
#' @param tau quantile to estimate
#' @param lambda penalty parameter
#' @param weights optional vector of weights
#' @param scale_x whether to scale the design matrix before estimation
#' @param method method argument to be passed to [quantreg::rq]
#' @param nfold number of folds to use when cross-validating
#' @param nlambda number of lambdas to search over when cross-validating
#' @param ... other arguments to pass to underlying fitting algorithm
rq.fit.post_lasso <- function(X, y, tau, lambda, weights,
                         scale_x = T, method = "agd", nfold = 5,
                         nlambda = 10, ...) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if(row(X) < nfold) {
    nfold = nrow(X)
  }

  intercept <- get_intercept(X)

  if(ncol(X) - (intercept > 0) <= 1) {
    stop("Model must have more than 1 covariate when using lasso.")
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weights)){
    message("No lambda provided, selecting based on ",nfold,"-fold cross-validation.")
    n = nrow(X)

    if(n != dim(as.matrix(weights))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weights

    # pre-multiplying the a matrix by a diagonal matrix of weights
    X <- diag(weights) %*% X
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
    message("No lambda provided, selecting based on ",nfold,"-fold cross-validation.")
    suppressWarnings({
      cv_fit <- lasso_cv_search(x = X[,-intercept], y = y , tau = tau,
                                method = method,
                                intercept = FALSE, nfolds = nfold,
                                foldid = NULL, nlambda = nlambda,
                                eps = 1e-04, init.lambda = 1,
                                parallel = T,
                                ...)
    })
    lambda = cv_fit$lambda.min
  }

  est <- fit_lasso(x = X[,-intercept],
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

  if(is.vector(new_X)) {
    new_X = matrix(new_X, ncol = 1)
  }

  est <- do.call(check_algorithm(method), args = list(X = new_X, y = y,
                                                      tau = tau,
                                                      method = method))
  coefficients[not_zero] <- est$coefficients

  residuals = y - X %*% coefficients

  if (!is.null(weights)){
    residuals <- residuals / weights
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = NA,
       ierr = 0,
       it = est$it,
       weights = weights,
       lambda = lambda)
}

#' Estimate a single quantile regression
#' @param X specification matrix for X variables
#' @param y outcome variable
#' @param tau quantile to regress
#' @param algorithm which algorithm to use
#' @param ... other arguments to be passed to the algorithm
#' @importFrom stats coefficients
#' @importFrom stats resid
#' @importFrom utils getFromNamespace
#' @importFrom methods existsFunction
fitQuantileRegression <- function(X, y, tau, algorithm = "rq.fit.sfn_start_val", ...) {
  qr_ns <- asNamespace("quantreg")
  qs_ns <- asNamespace("quantspace")

  if(algorithm %in% ls(qs_ns)) {
    f <- utils::getFromNamespace(algorithm, ns = "quantspace")
    do_matched_call(f, X = X, y = y, tau = tau, ...)
  } else if(algorithm %in% ls(qr_ns)) {

    # generically fit any of the quantreg algorithms
    args <- list(...)
    weights = args$weights

    if (!is.null(weights)){
      n = nrow(X)
      if(n != dim(as.matrix(weights))[1]){
        stop("Dimensions of design matrix and the weight vector not compatible")
      }

      # multiplying y by the weights
      y <- y * weights

      # pre-multiplying the a matrix by a diagonal matrix of weights
      X <- diag(weights) %*% X
    }

    if(!is.matrix(X)) {
      X <- as.matrix(X)
    }


    # wild hackery
    f <- utils::getFromNamespace(algorithm, "quantreg")

    fit <- do_matched_call(f, x = X, y = y, tau = tau, ...)

    list(coefficients = stats::coefficients(fit),
         residuals = stats::resid(fit),
         control = NA,
         ierr = 0,
         it = 0,
         weights = weights,
         lambda = list(...)$lambda)
  } else {
    f <- get(algorithm)
    X = as.matrix(X)

    do_matched_call(f, X = X, y = y, tau = tau, ...)
  }
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
#' @param weights vector of optional weights
#' @param ... other arguments to the function specified by the algorithm argument
#' @return List of estimated coefficients, warnings, iterations, and controls as in
#' standard quantile regression function
#' @export
regressResiduals = function(reg_spec_data,
                        ehat,
                        ind_hat,
                        tau,
                        trunc,
                        small,
                        weights,
                        algorithm = "rq.fit.sfn_start_val",
                        ...) {

  resids = ehat[ind_hat]

  if (!is.null(weights)){
    weights = as.matrix(weights)
  }

  if(trunc) {
    resids <- log(pmax(resids,small))
    spec_mat <- reg_spec_data$spec_matrix
  } else {
    weights = weights[resids > small]
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
    weights = as.vector(weights),
    ...)

  return(j_model)
}

#' return NA if argument is null
#' @param x value to check if null
na_if_null <- function(x) {
  if(is.null(x)){
      NA
  } else {
    x
  }
}


#' Computes coefficients for the quantile regression spacing method.
#' @param y Column of response variable.
#' @param X Regression specification matrix.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param small Minimum size of residuals for computational accuracy.
#' @param trunc Boolean value; if true, replace those dependent values less than small with small itself;
#' else, only use rows with residuals greater than small
#' @param algorithm The name of a function which will estimate a quantile regression.
#' Defaults to rq.fit.sfn_start_val. Must be a string, as it is passed to `do.call`
#' @param start_list Starting values for regression optimization.
#' @param weights vector of optional weights
#' @param control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`qs`]. This is set via the function [`qs_control`],
#' which returns a named list, with elements including:
#' * `trunc`: whether to truncate residual values below the argument "small"
#' * `small`: level of "small" values to guarentee numerical stability. If not specified, set dynamically based on the standard deviation of the outcome variable.
#' * `output_quantiles`: whether to save fitted quantiles as part of the function output
#' * `calc_avg_me`: whether to return average marginal effects as part of the fitted object
#' * `lambda`: the penalization factor to be passed to penalized regression algorithms
#' @param ... other parameters passed to the algorithm
#' @return
#' Returns a list of coefficients.
#' num_betas is an x by p matrix of estimated parameters for each supplied quantiles.
#' pseudo_r is a  1 by p matrix of psuedo R^2 values for each quantile estimate.
#' warnings is a 1 by p matrix of warnings produced by each quantile regression call.
#' iter: is a 1 by p matrix of iterations ran by each quantile regression call.
#' @export
#' @importFrom SparseM as.matrix
quantreg_spacing = function(
  y,
  X,
  var_names,
  alpha,
  jstar,
  algorithm = "rq.fit.sfn",
  weights = NULL,
  control = list(
    small = 1e-6,
    trunc = TRUE,
    start_list = NA,
    output_quantiles = FALSE,
    calc_avg_me = FALSE
  ),
  ...) {

  # disperse control parameters
  small = control$small
  trunc = control$trunc
  start_list = na_if_null(control$start_list)
  weights = weights
  output_quantiles = control$output_quantiles
  calc_avg_me = control$calc_avg_me
  lambda = control$lambda

  if(is.null(output_quantiles)) {
    output_quantiles <- TRUE
  }

  if(is.null(calc_avg_me)) {
    calc_avg_me <- FALSE
  }

  # number of observations that will default to an NA estimate of
  # a quantile.
  obs_thresh = 2

  if(nrow(X) < obs_thresh) {
    stop("Data contains fewer than ", obs_thresh, "observations.")
  }

  width = dim(X)[2]
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

  out_lambda = as.list(rep(NA, length(alpha)))

  # check to see if regression matrix is sparse. If not, then turn into CSR matrix
  if(!is(X, 'matrix.csr') & grepl("sfn", algorithm)) {
    X = denseMatrixToSparse(X)
  }

  if(missing(var_names) | is.null(var_names)) {
    var_names <- paste0("v", 1:ncol(X))
  }


  # Ensure matrix is not rank deficient
  reg_spec_starting_data <- list(spec_matrix = X, col_names = var_names)

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
      y = y,
      tau = tau,
      control = list( ),
      weights = weights,
      algorithm = algorithm,
      lambda = star_lambda,
      ...)
   } else {
     col_nums = getColNums(start_list, reg_spec_starting_data, alpha, jstar)
     sv = as.numeric(start_list[col_nums])
     star_model = fitQuantileRegression(
       X = reg_spec_starting_data$spec_matrix,
       y = y,
       tau = tau,
       control = list( ),
       sv = as.numeric(start_list[col_nums]),
       weights = weights,
       algorithm = algorithm,
       lambda = star_lambda,
       ...
      )
   }

  printWarnings(star_model)

  ehat0 = star_model$residuals

  #Calculate R^2
  V <- sum(rho(u = ehat0, tau = tau, weights = weights))
  V0 <- fitQuantileRegression(X = as.matrix.csr(rep(1, length(y))),
                             y = y,
                             tau = tau,
                             weights = weights,
                             algorithm = 'rq.fit.sfn')$residuals
  V0 <- sum(rho(u = V0, tau = tau,weights = weights))

  #set column names
  coef_df <- as.data.frame(t(star_model$coefficients))
  colnames(coef_df) <- reg_spec_starting_data$var_names
  # coef_df <- addMissingSpecColumns(
  #   coef_df,
  #   var_names)
  colnames(coef_df) <- paste(alpha[jstar], colnames(coef_df), sep="_")

  #log output for return
  pseudo_r[[jstar]] = (1 - V/V0)
  model[[jstar]] = coef_df
  warnings_log[[jstar]] = star_model$ierr
  iter_log[[jstar]] = star_model$it
  count_log[[jstar]] = dim(reg_spec_starting_data$spec_matrix)[1]
  out_lambda[[jstar]] = na_if_null(star_model$lambda)

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
        reg_spec_data <- list(spec_matrix = reg_spec_starting_data$spec_matrix[which(ehat > small),],
                                            col_names = reg_spec_starting_data$var_names)
      } else {
        # else, handle rank specification for typically-sized matrix
        reg_spec_data <- list(spec_matrix = reg_spec_starting_data$spec_matrix[ind_hat,],
                                            col_names = reg_spec_starting_data$var_names)
      }

      # run quantile regression
      coef <- NULL
      if(!any(is.na(start_list))){ # if user specified a start model, then use it as input
        col_nums = getColNums(start_list, reg_spec_data, alpha, j)
        sv = as.numeric(start_list[col_nums])
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    sv = sv, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                                    small = small,
                                    control = list( ),
                                    weights = weights[ind_hat],
                                    algorithm = algorithm, lambda = j_lambda,
                                    ...)
      } else{
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    ind_hat = ind_hat, tau = tau.t, trunc =  trunc,
                                    small = small,
                                    control = list( ),
                                    weights = weights[ind_hat],
                                    algorithm = algorithm, lambda = j_lambda,
                                    ...)
      }
      printWarnings(j_model)

      #Calculate R^2
      V <- sum(rho(u = j_model$residuals, tau = tau.t, weights = j_model$weights))
      V0 <- regressResiduals(
        reg_spec_data = list(
          spec_matrix = as.matrix.csr(rep(1, length(ind_hat)))
        ),
        ehat = ehat, ind_hat = ind_hat,
        tau = tau.t, trunc = trunc,
        small = small, weights = weights[ind_hat],
        algorithm = 'rq.fit.sfn')
      V0 <- sum(rho(u = V0$residuals, tau = tau.t,
                    weights = V0$weights))

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
      out_lambda[[j]] = na_if_null(j_model$lambda)

      # Update residuals
      ehat = ehat - exp(
        as.matrix(
          X %*%
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

      reg_spec_data <- list(spec_matrix = reg_spec_starting_data$spec_matrix[ind_hat,],
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
                                    control = list( ),
                                    algorithm = algorithm, lambda = j_lambda,
                                    weights = weights[ind_hat], ...)
      } else {
        j_model <- regressResiduals(reg_spec_data = reg_spec_data,
                                    ehat = -ehat,
                                    ind_hat = ind_hat,
                                    tau = tau.t, trunc = trunc, small = small,
                                    control = list( ),
                                    algorithm = algorithm,
                                    weights = weights[ind_hat],
                                    lambda = j_lambda,
                                    ...)
      }

      printWarnings(j_model)

      #Calculate pseudo-R^2
      V <- sum(rho(u = j_model$residuals, tau = tau.t, weights = j_model$weights))
      V0 <- regressResiduals(reg_spec_data = list('spec_matrix' = as.matrix.csr(rep(1, length(ind_hat)))),
                             ehat = -ehat, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                             small = small, weights = weights[ind_hat], algorithm = 'rq.fit.sfn')
      V0 <- sum(rho(u = V0$residuals, tau = tau.t, weights = V0$weights))

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
      out_lambda[[j]] = na_if_null(j_model$lambda)

        ehat = ehat + exp(as.matrix(
          X %*%
            unname(t(as.matrix(model[[j]])))))
    }
  }
  rv = list('coef' = do.call(cbind, model),
            'pseudo_r' = do.call(cbind, pseudo_r),
            'warnings' = do.call(cbind, warnings_log),
            'iter' = do.call(cbind, iter_log),
            'counts' = do.call(cbind, count_log),
            'lambda' = out_lambda)

  # calculate average marginal effects if user-specified
  if(calc_avg_me){
    qreg_coef = matrix(unlist(rv$coef), ncol = length(alpha))
    # calculate the average spacing (column-wise average)
    average_spacing = avg_spacing(X, qreg_coef, alpha,
                                    jstar, trim=.01)
    me = get_marginal_effects(qreg_coeffs = qreg_coef,
                              avg_spacings = average_spacing,
                              j_star = jstar,
                              calc_se = FALSE)
    rv$me = t(as.data.frame(as.vector(me$avgME)))
    colnames(rv$me) <- colnames(rv$coef)
    rownames(rv$me) <- NULL
  }

  # calculate quantiles induced by spacings if user-specified
  if(output_quantiles) {
    rv$quantiles = spacings_to_quantiles(spacingCoef = matrix(as.numeric(rv$coef),
                                                            ncol = p),
                                         X, jstar)
  }
  return(rv)
}

#' Get the data inside the s4 slot of this sparse matrix class
#' @param data some thing that could be a sparse matrix
#' @importFrom SparseM is.matrix.csr
get_underlying <- function(data) {
  if(SparseM::is.matrix.csr(data)) {
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

#' run function only with arguments
#' that match the arguments for f
#' @param f function
#' @param ... arguments to match
do_matched_call <- function(f, ...) {
  all_args = list(...)
  matched_args = subset(all_args,
                        names(all_args) %in%
                          names(formals(f)))
  do.call("f", args = matched_args)
}

