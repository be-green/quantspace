check <- function (x, tau = 0.5) {
  x * (tau - (x < 0))
}

#' Fit a quantile regression w/ a lasso penalty
#' @param x design matrix
#' @param y outcome variable
#' @param tau target quantile
#' @param lambda weight for penalization factor
#' @param weights optional observation weights
#' @param intercept whether or not to model the intercept
#' @param coef.cutoff what value for a coefficient is considered 0
#' @param method what underlying regression method to use for fitting
#' @param ... other arguments to pass to method
#' @importFrom SparseM as.matrix.csr
fit_lasso <- function (x, y, tau = 0.5, lambda = NULL, weights = NULL, intercept = TRUE,
          coef.cutoff = 1e-04, method = "sfn",
          ...) {
  if (is.null(dim(x))) {
    stop("x needs to be a matrix with more than 1 column")
  }
  p <- dim(x)[2]
  if (p == 1) {
    stop("x needs to be a matrix with more than 1 column")
  }
  n <- dim(x)[1]
  if (n != length(y)) {
    stop("length of y and rows of x do not match")
  }
  if (is.null(lambda) == TRUE | (length(lambda) != 1 & length(lambda) !=
                                 dim(x)[2])) {
    stop(paste("input of lambda must be of length 1 or",
               dim(x)[2]))
  }
  if (sum(lambda < 0) > 0) {
    stop(paste("lambda must be positive and we have a lambda of ",
               lambda, sep = ""))
  }

  lambda <- lambda * n
  if (length(lambda) == 1) {
    pen_x <- rbind(diag(rep(lambda, p)), diag(rep(-lambda,
                                                  p)))
  } else {
    pen_x <- rbind(diag(lambda), diag(-lambda))
    pen_x <- pen_x[rowSums(pen_x == 0) != dim(pen_x)[2],]
  }
  aug_n <- dim(pen_x)[1]
  aug_x <- rbind(x, pen_x)
  if (intercept) {
    aug_x <- cbind(c(rep(1, n), rep(0, aug_n)), aug_x)
  }
  aug_y <- c(y, rep(0, aug_n))

  if(!is.null(weights)) {
    orig_weights <- weights
    weights <- c(weights, rep(1, aug_n))

  }


  if(method == "sfn") {
    aug_x <- SparseM::as.matrix.csr(aug_x)
  }

  model = fitQuantileRegression(X = aug_x, y = aug_y,tau = tau,
                                algorithm = check_algorithm(method),
                                weights = weights, scale = 0,
                                ...)

  p_star <- p + intercept
  coefs <- coefficients(model)[1:p_star]
  return_val <- NULL
  return_val$coefficients <- coefs
  if (is.null(colnames(x))) {
    x_names <- paste("x", 1:p, sep = "")
  }
  else {
    x_names <- colnames(x)
  }
  if (intercept) {
    x_names <- c("intercept", x_names)
  }
  attributes(return_val$coefficients)$names <- x_names
  return_val$coefficients[abs(return_val$coefficients) < coef.cutoff] <- 0

  return_val$PenRho <- model$rho
  return_val$residuals <- model$residuals[1:n]

  if (is.null(weights)) {
    return_val$rho <- sum(sapply(return_val$residuals, check,
                                 tau))
  } else {
    return_val$rho <- sum(orig_weights * sapply(return_val$residuals,
                                                check, tau))
  }
  return_val$tau <- tau
  return_val$n <- n
  return_val$intercept <- intercept
  return_val
}

#' Randomly assign fold ids
#' @param n number of observations to assign
#' @param nfolds number of folds to assign
#' @details This is used internally for assignment of observations
#' for kfold cross-validation
randomly_assign <- function(n, nfolds) {
  group_size <- floor(n/nfolds)
  foldid <- rep(NA, n)
  obs <- 1:n
  for(i in 1:(nfolds - 1)) {
    group <- sample(obs, size = group_size, replace = F)
    foldid[group] <- i
    obs <- setdiff(obs, group)
  }
  foldid[obs] <- nfolds
  foldid
}

#' Search for optimal lambda via cross-validation
#' @param x design matrix for regression
#' @param y outcome variable for regression
#' @param tau target quantile
#' @param weights optional weights for regression
#' @param method method to be used for penalized quantile regression
#' (usually one of "sfn" or "br")
#' @param intercept Whether to model the intercept or not
#' @param nfolds number of folds to use for crossvalidation
#' @param eps smallest lambda used in search
#' @param nlambda number of lambdas to search over
#' @param foldid optional pre-specified fold identifier (for example,
#' if you want the folds to satisfy underlying data groupings)
#' @param init.lambda initial lambda for search
#' @param ... other parameters to pass on to fitting method
#' @param parallel whether to run cv scoring in parallel or not
#' @param coef.cutoff what cutoff to use for "0" coefficients
#' @param thresh threshhold for what counts as a "sparse enough"
#' solution for the top of the grid
#' @importFrom future.apply future_sapply
#' @importFrom future sequential
#' @importFrom future plan
#' @importFrom stats weighted.mean
#' @importFrom stats quantile
#' @details Searches for a range of lambdas by first expanding the
#' max lambda if necessary, then finding the smallest lambda
#' that sets all coefficients to zero. Then it computes kfold CV scores
#' along a grid of lambdas, returning the scores and the smallest lambda.
lasso_cv_search <- function (x, y, tau = 0.5,
                             weights = NULL, method = "sfn",
                             intercept = TRUE, nfolds = 10,
                             foldid = NULL, nlambda = 100,
                             eps = 1e-04, init.lambda = 2,
                             parallel = T, coef.cutoff = 1e-5,
                             thresh = 0.01,
                             ...) {
  p <- dim(x)[2]

  p_range <- 1:p + intercept
  n <- dim(x)[1]
  # penalty function is lasso penalty
  pen_func <- function(x, lambda) lambda * abs(x)

  lambda_star <- init.lambda
  searching <- TRUE
  while (searching) {

    init_fit <- fit_lasso(x, y, tau, lambda = lambda_star,
                          method = method,
                          weights = weights,
                          intercept = intercept, ...)

    if (mean(abs(init_fit$coefficients[p_range])) <= thresh) {
      searching <- FALSE
    } else {
      lambda_star <- lambda_star * 1.5
    }
  }
  lambda_min <- eps * lambda_star
  lambda <- exp(seq(log(max(lambda_min)), log(max(lambda_star)),
                    length.out = nlambda))
  models <- list()
  fit_models <- TRUE

  # bracket search until it finds set of lambdas that don't always zero out
  # coefficients

  # start in the middle
  lam_pos = floor(length(lambda)/2)
  max_lambda_pos <- length(lambda)
  min_lambda_pos <- 1
  while (fit_models & lambda[max_lambda_pos] > 2) {
    if (fit_models) {
      models[[lam_pos]] <- fit_lasso(x = x,  y = y, tau = tau,
                                     lambda = lambda[lam_pos],
                                     weights =  weights,
                                     intercept = intercept,
                                     method = method,
                                     ...
                                     )
    }
    if (mean(abs(coefficients(models[[lam_pos]])[p_range])) > thresh) {
      if(lam_pos > min_lambda_pos) {
        min_lambda_pos = lam_pos
      }
      lam_pos = lam_pos + ceiling((max_lambda_pos - lam_pos)/2)
    } else {
      if(lam_pos == max_lambda_pos) {
        end_pos = lam_pos
        break
      }
      if(lam_pos < max_lambda_pos) {
        max_lambda_pos = lam_pos
      }
      lam_pos = lam_pos - floor((lam_pos - min_lambda_pos)/2)
    }
  }

  lambda = lambda[1:lam_pos]

  if (is.null(foldid)) {
    foldid <- randomly_assign(n, nfolds)
  }

  # generate all combinations of lambda and folds
  lambda_folds <- expand.grid(lambda, 1:nfolds)
  colnames(lambda_folds) <- c("lambda", "foldid")

  fold_list <- list()
  for(i in 1:nrow(lambda_folds)) {
    fold_list[[i]] <- lambda_folds[i,]
  }

  if(parallel) {
    setCores(getCores())
  } else {
    future::plan(future::sequential)
  }

  cv_results <-
    future.apply::future_sapply(fold_list,
                                FUN = function(l, x, y, foldid) {

    lambda = l[[1]]
    i = l[[2]]
    train_x <- x[foldid != i, ]
    train_y <- y[foldid != i]
    test_x <- x[foldid == i, , drop = FALSE]
    test_y <- y[foldid == i]

    if (is.null(weights)) {
      train_weights <- test_weights <- NULL
    }   else {
      train_weights <- weights[foldid != i]
      test_weights <- weights[foldid == i]
    }

    model <- fit_lasso(x = train_x, y = train_y, tau = tau,
                       weights = train_weights,
                       intercept = intercept,
                       method = method, lambda = lambda)

    model_coef <- model$coefficients

    # if test_x matrix doesn't have intercept, add it
    if(ncol(test_x) == (length(model_coef) - 1)) {
      test_x = cbind(1, test_x)
    }

    y_hat <- test_x %*% model_coef

    model_fit <- check(test_y - y_hat)

    if(is.null(test_weights)) {
      mean(model_fit)
    } else {
      stats::weighted.mean(model_fit, test_weights)
    }
  }, x = x, y = y, foldid = foldid)


  cv_results <- sapply(split(cv_results, factor(lambda_folds$lambda)), mean)
  lambda.min <- lambda[which.min(cv_results)]
  return_val <- list()
  return_val$models <- models
  return_val$cv <- data.frame(lambda = lambda, cve = cv_results)
  rownames(return_val$cv) <- NULL
  colnames(return_val$cv)[2] <- "cv_score"
  return_val$lambda.min <- lambda.min

  return_val
}
