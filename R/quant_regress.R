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
#' @export
rq.fit.sfn_start_val <- function(X,y,tau=.5,
                                 rhs = (1-tau)*c(t(a) %*% rep(1,length(y))),
                                 control,
                                 sv,
                                 weight_vec = NULL) {
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
  ctrl$small <- ctrl$small/100
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

  residuals <- -y - a %*% coefficients
  if (!is.null(weight_vec)){
    residuals <- residuals / weight_vec
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = ctrl,
       ierr = ierr,
       it = fit$maxiter,
       weight_vec = weight_vec)
}

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
  reordered_coefs <- rep(NA, length(coefficients))
  reordered_coefs[-intercept] <- coefficients[-1]
  reordered_coefs[intercept] <- coefficients[1]
  reordered_coefs
}

#' Quantile Regression w/ Lasso Penalty
#' @param X Design matrix, X
#' @param y outcome variable, y
#' @param tau quantile to estimate
#' @param lambda penalty parameter
#' @param weight_vec optional vector of weights
#' @param scale_x whether to scale the design matrix before estimation
#' @param ... other arguments to pass to rqPen::rq.lasso.fit
#' @importFrom rqPen rq.lasso.fit
rq.fit.lasso <- function(X, y, tau, lambda, weight_vec, scale_x = T, ...) {

  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }
  intercept <- get_intercept(X)

  # get the intercept, and then scale everything except the intercept
  if (scale_x) {
    unscaled_X <- X
    X <- scale_for_lasso(X, intercept)
    mu_x <- attr(X, "scaled:center")
    sigma_x <- attr(X, "scaled:scale")
    attr(X, "scaled:center") <- NULL
    attr(X, "scaled:scale") <- NULL
  }

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){
    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }

    # multiplying y by the weights
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    X <- diag(weight_vec) %*% X
  }

  est <- rqPen::rq.lasso.fit(x = X[,-intercept],
                             y = y,
                             tau = tau,
                             lambda = lambda,
                             intercept = T,
                             scalex = F, ...)

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
       weight_vec = weight_vec)
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
    weight_vec = as.matrix(weight_vec[ind_hat])
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
    weight_vec,
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
#' @param ... other parameters passed to the algorithm
#' @import SparseM
#' @return
#' Returns a list of coefficients.
#' num_betas is an x by p matrix of estimated parameters for each supplied quantiles.
#' pseudo_r is a  1 by p matrix of psuedo R^2 values for each quantile estimate.
#' warnings is a 1 by p matrix of warnings produced by each quantile regression call.
#' iter: is a 1 by p matrix of iterations ran by each quantile regression call.
#' @export
quantRegSpacing = function(
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
  ...) {

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

  # Calculate initial fit
  if(any(is.na(start_list))){ # if user supplied starting values, then use them
    star_model = fitQuantileRegression(
      X = reg_spec_starting_data$spec_matrix,
      y = dep_col,
      tau = tau,
      control = list(tmpmax = tmpmax),
      weight_vec = weight_vec,
      algorithm = algorithm,
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

  # Estimate upper quantiles sequentially
  ehat = ehat0
  for (j in (jstar+1):p) {
    ind_hat = which(ehat > 0)

    if(length(ind_hat) == 0) {
      # Update residuals
      coef = rep(NA, length(star_model$coefficients))
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


      #run quantile regression
      coef <- NULL
      if(!any(is.na(start_list))){ # if user specified a start model, then use it as input
        col_nums = getColNums(start_list, reg_spec_data, alpha, j)
        sv = as.numeric(start_list[col_nums])
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    sv = sv, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                                    small = small, control = list(tmpmax = tmpmax),
                                    weight_vec = weight_vec,
                                    algorithm = algorithm, ...)
      } else{
        j_model <- regressResiduals(reg_spec_data = reg_spec_data, ehat = ehat,
                                    ind_hat = ind_hat, tau = tau.t, trunc =  trunc,
                                    small = small, control = list(tmpmax = tmpmax),
                                    weight_vec = weight_vec,
                                    algorithm = algorithm,
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
        small = small, weight_vec = weight_vec,
        algorithm = 'rq.fit.sfn_start_val')
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

    # Ensure the cut of the starting data that we take for
    # current spacing is not rank-deficient
    # if not truncating, then ensure exact rows of the regression matrix is included
    # else, handle rank specification for typically-sized matrix
    if(!trunc) {
      ind_hat = which(-ehat > small)
    } else {
      ind_hat = which(ehat < 0)
    }

    if(length(ind_hat) == 0) {
      # Update residuals
      coef = rep(NA, length(star_model$coefficients))
      coef_df <- as.data.frame(t(coef))

      #get column names
      colnames(coef_df) <- reg_spec_data$var_names
      coef_df <- addMissingSpecColumns(
        coef_df,
        var_names)
      colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

      #log results
      model[[j]] = coef_df
      pseudo_r[[j]] = NA
      warnings_log[[j]] = j_model$ierr
      iter_log[[j]] = j_model$it
      count_log[[j]] = dim(reg_spec_data$spec_matrix)[1]

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
                                    algorithm = algorithm,
                                    weight_vec = weight_vec, ...)
      } else {
        j_model <- regressResiduals(reg_spec_data = reg_spec_data,
                                    ehat = -ehat,
                                    ind_hat = ind_hat,
                                    tau = tau.t, trunc = trunc, small = small,
                                    control = list(tmpmax = tmpmax),
                                    algorithm = algorithm,
                                    weight_vec = weight_vec, ...)
      }

      printWarnings(j_model)

      #Calculate pseudo-R^2
      V <- sum(rho(u = j_model$residuals, tau = tau.t, weight_vec = j_model$weight_vec))
      V0 <- regressResiduals(reg_spec_data = list('spec_matrix' = as.matrix.csr(rep(1, length(ind_hat)))),
                             ehat = -ehat, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                             small = small, weight_vec = weight_vec, algorithm = 'rq.fit.sfn_start_val')
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

        ehat = ehat + exp(as.matrix(
          data %*%
            unname(t(as.matrix(model[[j]])))))
    }
  }
  rv = list('coef' = do.call(cbind, model),
            'pseudo_r' = do.call(cbind, pseudo_r),
            'warnings' = do.call(cbind, warnings_log),
            'iter' = do.call(cbind, iter_log),
            'counts' = do.call(cbind, count_log))

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
    rv$quantiles = spacingsToQuantiles(spacingCoef = matrix(as.numeric(rv$coef),
                                                            ncol = p),
                                       data, jstar)
  }
  return(rv)
}

#' Compute quantiles given parameter coefficients and data
#' @param spacingCoef J by p matrix; row is number of variables, p is number of quantiles
#' @param data independent variables
#' @param jstar index of median quantiles
#' @return N by p matrix of quantiles
#' @export
spacingsToQuantiles <- function(spacingCoef, data, jstar) {

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
