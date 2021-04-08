#' Function which gets resampled rows
#' @param rows rows to resample
#' @param sampling_method method for resampling
#' @param ... other arguments passed on to the method
getRows <- function(rows, sampling_method, ...) {
  do.call(sampling_method, args = list(rows = rows, ...))
}


#' Function that doesn't reorder rows (for weighted bootstrap)
#' @param rows rows to potentially resample
#' @param ... other arguments, ignored
leaveRows <- function(rows, ...) {
  rows
}

#' Function that fully bootstraps rows
#' @param rows rows to potentially resample
#' @param ... other arguments, ignored
bootstrapRows <- function(rows, ...) {
  sample(rows, replace = T)
}

#' Function that fully bootstraps rows
#' @param rows rows to potentially resample
#' @param subsample_percent percentage of data to use in subsampling
#' @param ... other arguments, ignored
subsampleRows <- function(rows, subsample_percent, ...) {
  sample(rows,floor(subsample_percent*length(rows)), replace = F)
}

#' Function that returns the correct weights for weighted bootstrap
#' @param weights vector of existing regression weights
#' @param draw_weights whether draw random exponential weights
#' @param n number of weights to draw
getWeights <- function(weights, draw_weights, n) {
  if(draw_weights) {
    new_weights <- pmax(stats::rexp(n),1e-3)
    if(is.null(weights)) {
      new_weights
    } else {
      new_weights * weights
    }
  } else {
    weights
  }
}

#' @rdname standard_errors
#' @param subsample_percent A number between 0 and one, specifying the percent of the data to subsample for standard error calculations
#' @param draw_weights Whether to use random exponential weights for bootstrap, either TRUE or FALSE
#' @param sampling_method One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' @param ... Other arguments, ignored for now
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and bootstrapRows
#' samples with replacement.
#' @export
resample_qs <- function(X,
                        y,
                        weights,
                        sampling_method,
                        alpha,
                        jstar,
                        control,
                        algorithm,
                        draw_weights,
                        var_names,
                        subsample_percent,
                        ...) {


  if(!grepl("Rows$", sampling_method)) {
    sampling_method <- paste0(sampling_method, "Rows")
  }

  rows <- getRows(1:nrow(X), sampling_method = sampling_method, subsample_percent = subsample_percent)
  weights <- getWeights(weights, draw_weights, length(rows))

  quantreg_spacing(
    y = y[rows],
    X = X[rows,],
    var_names = var_names,
    alpha = alpha,
    jstar = jstar,
    control = control,
    weights = weights,
    algorithm = algorithm,
    ...)
}

#' @rdname standard_errors
#' @export
weighted_bootstrap <- function(X,
                               y,
                               weights,
                               sampling_method,
                               alpha,
                               jstar,
                               control,
                               algorithm,
                               draw_weights,
                               var_names,
                               subsample_percent,
                               ...) {
  resample_qs(X = X, y = y, weights = weights,
           sampling_method = "leaveRows", draw_weights = T,
           alpha = alpha, jstar = jstar, control = control,
           algorithm = algorithm, subsample_percent = 1,
           var_names = var_names, ...)
}

#' @rdname standard_errors
#' @export
bootstrap <- function(X,
                      y,
                      weights,
                      sampling_method,
                      alpha,
                      jstar,
                      control,
                      algorithm,
                      draw_weights,
                      var_names,
                      subsample_percent,
                      ...) {
  resample_qs(X = X, y = y, weights = weights,
              sampling_method = "bootstrapRows", draw_weights = F,
              alpha = alpha, jstar = jstar, control = control,
              algorithm = algorithm, subsample_percent = 1,
              var_names = var_names, ...)
}

#' @rdname standard_errors
subsample <- function(X,
                      y,
                      weights,
                      sampling_method,
                      alpha,
                      jstar,
                      control,
                      algorithm,
                      draw_weights,
                      var_names,
                      subsample_percent,
                      ...) {
  resample_qs(X = X, y = y, weights = weights,
              sampling_method = "subsampleRows", draw_weights = F,
              alpha = alpha, jstar = jstar, control = control,
              algorithm = algorithm, subsample_percent = subsample_percent,
              var_names = var_names, ...)
}

#' check if all values of two vectors match
#' @param x first vector
#' @param y second vector
#' @importFrom assertthat assert_that
all_match <- function(x, y) {
  all(x == y)
}

#' Check which rows of a matrix match a whole vector
#' @param X first matrix
#' @param Y vector to match
#' @importFrom assertthat assert_that
columns_match_vector <- function(X, Y) {
  assertthat::assert_that(ncol(X) == length(Y))

  matched_cols <- vector()

  for(i in 1:nrow(X)) {
    matched_cols[i] <- all_match(X[i, ], Y)
  }

  matched_cols

}

#' Get clusters for subsampling given a formula
#' @param cluster_matrix matrix of cluster column values
get_strata <- function(cluster_matrix) {

  if(is.vector(cluster_matrix)) {
    cluster_matrix <- as.matrix(cluster_matrix, ncol = 1)
  }

  if(is.data.frame(cluster_matrix)) {
    cluster_matrix <- as.matrix(cluster_matrix)
  }

  levels <- unique(cluster_matrix)
  level_vec <- rep(NA, length(levels))

  for(i in 1:nrow(levels)) {
    level_vec[which(columns_match_vector(cluster_matrix, levels[i,]))] <- i
  }

  level_vec
}

#' Computes standard errors for the quantile regression spacing method using
#' subsampling.
#' @param X Regression specification matrix.
#' @param y Column of response variable.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param parallel whether to run in parallel or not
#' @param control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`qs`]. This is a named list, with possible elements including:
#' * `trunc`: whether to truncate residual values below the argument "small"
#' * `small`: level of "small" values to guarentee numerical stability. If not specified, set dynamically based on the standard deviation of the outcome variable.
#' * `lambda`: For penalized regression, you can specify a level of lambda which will weight the penalty. If not set, will be determined based on 10-fold cross-validation.
#' * `output_quantiles`: whether to save fitted quantiles as part of the function output
#' * `calc_avg_me`: whether to return average marginal effects as part of the fitted object
#' * Any other arguments passed to specific downstream quantile regression algorithms (e.g. rq.fit).
#' @param cluster_matrix Matrix of cluster variables, as returned by a model formula
#' @param std_err_control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`standard_errors`]. Possible arguments include:
#' * `se_method`: Method to use for standard errors, either "weighted_bootstrap",
#' "subsample", "bootstrap" or "custom" along with a specified subsampling method and
#' subsample percent. If specifying "custom", must also specify `subsampling_percent` and
#' `draw_weights`. If you specify "subsample", subsampling percent defaults to 0.2, but can be
#' changed. See details for details.
#' * `num_bs`: Number of bootstrap iterations to use, defaults to 100.
#' * `subsample_percent`: A number between 0 and one, specifying the percent of the data to subsample for standard error calculations
#' * `draw_weights`: Whether to use random exponential weights for bootstrap, either TRUE or FALSE
#' * `sampling_method` One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and boostrapRows
#' samples with replacement.`
#' @param weights vector of same length and order as dependent column, to be used as weights for estimation
#'              (note, if draw weights is set to TRUE, this variable will be the element-wise product
#'              of itself and a random vector of weights)
#' @param algorithm function which is actually used to fit each quantile regression
#' @param seed Seed to be used when generating RNG
#' @param ... other arguments passed to quantile fitting function
#' @importFrom methods is
#' @importFrom stats rexp
#' @importFrom stats cov
#' @importFrom purrr pmap
#' @importFrom stats dnorm
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future.apply future_lapply
#' @importFrom purrr quietly
#' @export
standard_errors <- function(y,
                            X,
                            cluster_matrix,
                            algorithm,
                            control = qs_control(),
                            std_err_control = std_err_control(),
                            var_names,
                            alpha,
                            jstar,
                            parallel = F,
                            weights = NULL,
                            seed = NULL,
                            ...) {

  # if they set parallel to false, we still
  # want to restore the old plan afterwards
  if(!parallel) {
    old_plan = future::plan()
    future::plan(future::sequential)
  }

  if(is.null(seed)) {
    seed = T
  }

  if(is.null(cluster_matrix)) {
    cluster_index <- rep(1, nrow(X))
  } else {
    cluster_index <- get_strata(cluster_matrix)
  }

  if(std_err_control$se_method %in% c("bootstrap", "weighted_bootstrap")) {
    subsample_percent = 1
  } else {
    subsample_percent = std_err_control$subsample_percent
  }

  num_clusters <- length(unique(cluster_index))
  num_bs <- std_err_control$num_bs

  fits <- future.apply::future_lapply(1:num_bs,
                                      future.seed = seed,
                            purrr::quietly(
                            function(x, X_mat,
                                     y,
                                     cluster_index,
                                     num_clusters,
                                     algorithm,
                                     se_method,
                                     var_names,
                                     alpha,
                                     jstar,
                                     control,
                                     weights,
                                     subsample_percent,
                                     draw_weights,
                                     sampling_method = sampling_method,
                                     ...) {
                              clustered_fits <- list()
                              for(j in 1:num_clusters) {
                                clustered_fits[[j]] <-
                                  do.call(se_method,
                                        args = list(X = X_mat[which(cluster_index == j),],
                                                    y = y[which(cluster_index == j)],
                                                    weights = weights[which(cluster_index == j)],
                                                    alpha = alpha, jstar = jstar, control = control,
                                                    algorithm = algorithm,
                                                    subsample_percent = subsample_percent,
                                                    var_names = var_names,
                                                    draw_weights = draw_weights,
                                                    sampling_method = sampling_method,
                                                    ...))
                              }
                              purrr::pmap(clustered_fits, rbind)
  }), X_mat = X, se_method = std_err_control$se_method, cluster_index = cluster_index,
  y = y, weights = weights, alpha = alpha,
  jstar = jstar, var_names = var_names, control = control,
  algorithm = algorithm,
  subsample_percent = subsample_percent, sampling_method = std_err_control$sampling_method,
  num_clusters = num_clusters, draw_weights = std_err_control$draw_weights,  ...)

  warns <- Filter(function(x) length(x) > 0, lapply(fits, function(x) x$warnings))

  if(length(warns) > 0) {
    if(length(unique(unlist(warns))) > 10) {
      warning("\n",length(warns), " out of ", num_bs, " samples produced warnings. Showing first 10:\n\t",
              paste0(head(unique(unlist(warns)), 10), collapse = "\n\t"))
    } else {
      warning("\n",length(warns), " out of ", num_bs, " samples produced warnings:\n\n",
              paste0(unique(unlist(warns)), collapse = "\n"))
    }


  }

  fit <- purrr::pmap(lapply(fits, function(x) x$result), rbind)
  quant_cov_mat <- stats::cov(fit$coef, use = "pairwise.complete.obs")
  quant_cov_mat <- quant_cov_mat * subsample_percent

  if(!parallel) {
    future::plan(old_plan)
  }


  return(list('quant_cov' = quant_cov_mat,
              'warnings' = fit$warnings,
              'iter' = fit$iter,
              'counts' = fit$counts,
              'coef_boot' = fit$coef,
              'subsample_percent' = subsample_percent))
}
#' Function which gets resampled rows
#' @param rows rows to resample
#' @param sampling_method method for resampling
#' @param ... other arguments passed on to the method
getRows <- function(rows, sampling_method, ...) {
  do.call(sampling_method, args = list(rows = rows, ...))
}


#' Function that doesn't reorder rows (for weighted bootstrap)
#' @param rows rows to potentially resample
#' @param ... other arguments, ignored
leaveRows <- function(rows, ...) {
  rows
}

#' Function that fully bootstraps rows
#' @param rows rows to potentially resample
#' @param ... other arguments, ignored
bootstrapRows <- function(rows, ...) {
  sample(rows, replace = T)
}

#' Function that fully bootstraps rows
#' @param rows rows to potentially resample
#' @param subsample_percent percentage of data to use in subsampling
#' @param ... other arguments, ignored
subsampleRows <- function(rows, subsample_percent, ...) {
  sample(rows,floor(subsample_percent*length(rows)), replace = F)
}

#' Function that returns the correct weights for weighted bootstrap
#' @param weights vector of existing regression weights
#' @param draw_weights whether draw random exponential weights
#' @param n number of weights to draw
getWeights <- function(weights, draw_weights, n) {
  if(draw_weights) {
    new_weights <- pmax(stats::rexp(n),1e-3)
    if(is.null(weights)) {
      new_weights
    } else {
      new_weights * weights
    }
  } else {
    weights
  }
}

#' @rdname standard_errors
#' @export
resample_qs <- function(X,
                        y,
                        weights,
                        sampling_method,
                        alpha,
                        jstar,
                        control,
                        algorithm,
                        draw_weights,
                        var_names,
                        subsample_percent,
                        ...) {


  if(!grepl("Rows$", sampling_method)) {
    sampling_method <- paste0(sampling_method, "Rows")
  }

  rows <- getRows(1:nrow(X), sampling_method = sampling_method, subsample_percent = subsample_percent)
  weights <- getWeights(weights, draw_weights, length(rows))

  quantreg_spacing(
    y = y[rows],
    X = X[rows,],
    var_names = var_names,
    alpha = alpha,
    jstar = jstar,
    control = control,
    weights = weights,
    algorithm = algorithm,
    ...)
}

#' @rdname standard_errors
#' @export
weighted_bootstrap <- function(X,
                               y,
                               weights,
                               sampling_method,
                               alpha,
                               jstar,
                               control,
                               algorithm,
                               draw_weights,
                               var_names,
                               subsample_percent,
                               ...) {
  resample_qs(X = X, y = y, weights = weights,
           sampling_method = "leaveRows", draw_weights = T,
           alpha = alpha, jstar = jstar, control = control,
           algorithm = algorithm, subsample_percent = 1,
           var_names = var_names, ...)
}

#' @rdname standard_errors
#' @export
bootstrap <- function(X,
                      y,
                      weights,
                      sampling_method,
                      alpha,
                      jstar,
                      control,
                      algorithm,
                      draw_weights,
                      var_names,
                      subsample_percent,
                      ...) {
  resample_qs(X = X, y = y, weights = weights,
              sampling_method = "bootstrapRows", draw_weights = F,
              alpha = alpha, jstar = jstar, control = control,
              algorithm = algorithm, subsample_percent = 1,
              var_names = var_names, ...)
}

#' @rdname standard_errors
subsample <- function(X,
                      y,
                      weights,
                      sampling_method,
                      alpha,
                      jstar,
                      control,
                      algorithm,
                      draw_weights,
                      var_names,
                      subsample_percent,
                      ...) {
  resample_qs(X = X, y = y, weights = weights,
              sampling_method = "subsampleRows", draw_weights = F,
              alpha = alpha, jstar = jstar, control = control,
              algorithm = algorithm, subsample_percent = subsample_percent,
              var_names = var_names, ...)
}

#' check if all values of two vectors match
#' @param x first vector
#' @param y second vector
#' @importFrom assertthat assert_that
all_match <- function(x, y) {
  all(x == y)
}

#' Check which rows of a matrix match a whole vector
#' @param X first matrix
#' @param Y vector to match
#' @importFrom assertthat assert_that
columns_match_vector <- function(X, Y) {
  assertthat::assert_that(ncol(X) == length(Y))

  matched_cols <- vector()

  for(i in 1:nrow(X)) {
    matched_cols[i] <- all_match(X[i, ], Y)
  }

  matched_cols

}

#' Get clusters for subsampling given a formula
#' @param cluster_matrix matrix of cluster column values
get_strata <- function(cluster_matrix) {

  if(is.vector(cluster_matrix)) {
    cluster_matrix <- as.matrix(cluster_matrix, ncol = 1)
  }

  if(is.data.frame(cluster_matrix)) {
    cluster_matrix <- as.matrix(cluster_matrix)
  }

  levels <- unique(cluster_matrix)
  level_vec <- rep(NA, length(levels))

  for(i in 1:nrow(levels)) {
    level_vec[which(columns_match_vector(cluster_matrix, levels[i,]))] <- i
  }

  level_vec
}

single_bootstrap_iter <- function(x,
                                  X_mat,
                                  y,
                                  cluster_index,
                                  num_clusters,
                                  algorithm,
                                  se_method,
                                  var_names,
                                  alpha,
                                  jstar,
                                  control,
                                  weights,
                                  subsample_percent,
                                  draw_weights,
                                  sampling_method,
                                  ...) {
  clustered_fits <- list()
  for(j in 1:num_clusters) {
    clustered_fits[[j]] <-
      do.call(se_method,
              args = list(X = X_mat[which(cluster_index == j),],
                          y = y[which(cluster_index == j)],
                          weights = weights[which(cluster_index == j)],
                          alpha = alpha, jstar = jstar, control = control,
                          algorithm = algorithm,
                          subsample_percent = subsample_percent,
                          var_names = var_names,
                          draw_weights = draw_weights,
                          sampling_method = sampling_method,
                          ...))
  }
  purrr::pmap(clustered_fits, rbind)
}

#' Computes standard errors for the quantile regression spacing method using
#' subsampling.
#' @param X Regression specification matrix.
#' @param y Column of response variable.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param parallel whether to run in parallel or not
#' @param control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`qs`]. This is a named list, with possible elements including:
#' * `trunc`: whether to truncate residual values below the argument "small"
#' * `small`: level of "small" values to guarentee numerical stability. If not specified, set dynamically based on the standard deviation of the outcome variable.
#' * `lambda`: For penalized regression, you can specify a level of lambda which will weight the penalty. If not set, will be determined based on 10-fold cross-validation.
#' * `output_quantiles`: whether to save fitted quantiles as part of the function output
#' * `calc_avg_me`: whether to return average marginal effects as part of the fitted object
#' * Any other arguments passed to specific downstream quantile regression algorithms (e.g. rq.fit).
#' @param cluster_matrix Matrix of cluster variables, as returned by a model formula
#' @param std_err_control control parameters to pass to the control arguments of [`quantreg_spacing`],
#' the lower-level function called by [`standard_errors`]. Possible arguments include:
#' * `se_method`: Method to use for standard errors, either "weighted_bootstrap",
#' "subsample", "bootstrap" or "custom" along with a specified subsampling method and
#' subsample percent. If specifying "custom", must also specify `subsampling_percent` and
#' `draw_weights`. If you specify "subsample", subsampling percent defaults to 0.2, but can be
#' changed. See details for details.
#' * `num_bs`: Number of bootstrap iterations to use, defaults to 100.
#' * `subsample_percent`: A number between 0 and one, specifying the percent of the data to subsample for standard error calculations
#' * `draw_weights`: Whether to use random exponential weights for bootstrap, either TRUE or FALSE
#' * `sampling_method` One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and boostrapRows
#' samples with replacement.`
#' @param weights vector of same length and order as dependent column, to be used as weights for estimation
#'              (note, if draw weights is set to TRUE, this variable will be the element-wise product
#'              of itself and a random vector of weights)
#' @param algorithm function which is actually used to fit each quantile regression
#' @param seed Seed to be used when generating RNG
#' @param ... other arguments passed to quantile fitting function
#' @importFrom methods is
#' @importFrom stats rexp
#' @importFrom stats cov
#' @importFrom purrr pmap
#' @importFrom stats dnorm
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future.apply future_lapply
#' @importFrom purrr quietly
#' @export
standard_errors <- function(y,
                            X,
                            cluster_matrix,
                            algorithm,
                            control = qs_control(),
                            std_err_control = std_err_control(),
                            var_names,
                            alpha,
                            jstar,
                            parallel = F,
                            weights = NULL,
                            seed = NULL,
                            ...) {


  num_bs <- std_err_control$num_bs


  if(is.null(seed)) {
    seed = T
  }

  if(is.null(cluster_matrix)) {
    cluster_index <- rep(1, nrow(X))
  } else {
    cluster_index <- get_strata(cluster_matrix)
  }

  if(std_err_control$se_method %in% c("bootstrap", "weighted_bootstrap")) {
    subsample_percent = 1
  } else {
    subsample_percent = std_err_control$subsample_percent
  }

  num_clusters <- length(unique(cluster_index))


  parallel_thresh = std_err_control$parallel_thresh

  # if they set parallel to false, we still
  # want to restore the old plan afterwards
  if(num_bs * nrow(X) * num_clusters < parallel_thresh) {
    parallel = FALSE
  }

  if(!parallel) {
    old_plan = future::plan()
    future::plan(future::sequential)
  }

  fits <- future.apply::future_lapply(1:num_bs,
                                      future.seed = seed,
                            purrr::quietly(single_bootstrap_iter),
                            X_mat = X, se_method = std_err_control$se_method, cluster_index = cluster_index,
  y = y, weights = weights, alpha = alpha,
  jstar = jstar, var_names = var_names, control = control,
  algorithm = algorithm,
  subsample_percent = subsample_percent, sampling_method = std_err_control$sampling_method,
  num_clusters = num_clusters, draw_weights = std_err_control$draw_weights,  ...)

  warns <- Filter(function(x) length(x) > 0, lapply(fits, function(x) x$warnings))

  if(length(warns) > 0) {
    if(length(unique(unlist(warns))) > 10) {
      warning("\n",length(warns), " out of ", num_bs, " samples produced warnings. Showing first 10:\n\t",
              paste0(head(unique(unlist(warns)), 10), collapse = "\n\t"))
    } else {
      warning("\n",length(warns), " out of ", num_bs, " samples produced warnings:\n\n",
              paste0(unique(unlist(warns)), collapse = "\n"))
    }


  }

  fit <- purrr::pmap(lapply(fits, function(x) x$result), rbind)
  quant_cov_mat <- stats::cov(fit$coef, use = "pairwise.complete.obs")
  quant_cov_mat <- quant_cov_mat * subsample_percent

  if(!parallel) {
    future::plan(old_plan)
  }


  return(list('quant_cov' = quant_cov_mat,
              'warnings' = fit$warnings,
              'iter' = fit$iter,
              'counts' = fit$counts,
              'coef_boot' = fit$coef,
              'subsample_percent' = subsample_percent))
}
