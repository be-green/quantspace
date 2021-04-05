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
#' @param weight_vec vector of existing regression weights
#' @param draw_weights whether draw random exponential weights
#' @param n number of weights to draw
getWeights <- function(weight_vec, draw_weights, n) {
  if(draw_weights) {
    new_weights <- pmax(stats::rexp(n),1e-3)
    if(is.null(weight_vec)) {
      new_weights
    } else {
      new_weights * weight_vec
    }
  } else {
    weight_vec
  }
}

#' @rdname standard_errors
#' @export
resample_qs <- function(data,
                       dep_col,
                       weight_vec,
                       sampling_method,
                       draw_weights,
                       alpha,
                       jstar,
                       small,
                       trunc,
                       algorithm,
                       subsample_percent,
                       var_names,
                       ...) {


  if(!grepl("Rows$", sampling_method)) {
    sampling_method <- paste0(sampling_method, "Rows")
  }

  rows <- getRows(1:nrow(data), sampling_method = sampling_method, subsample_percent = subsample_percent)
  weights <- getWeights(weight_vec, draw_weights, length(rows))

  quantRegSpacing(
    data = data[rows,],
    dep_col = dep_col[rows],
    var_names = var_names,
    alpha = alpha,
    jstar = jstar,
    small = small,
    trunc = trunc,
    start_list = NA,
    weight_vec = weights,
    algorithm = algorithm,
    ...)
}

#' @rdname standard_errors
#' @export
weighted_bootstrap <- function(data,
                               dep_col,
                               weight_vec,
                               sampling_method,
                               alpha,
                               jstar,
                               small,
                               trunc,
                               algorithm,
                               draw_weights,
                               var_names,
                               subsample_percent,
                               ...) {
  resample_qs(data, dep_col, weight_vec,
           sampling_method = "leaveRows", draw_weights = T,
           alpha, jstar, small, trunc, algorithm, subsample_percent = 1,
           var_names = var_names, ...)
}

#' @rdname standard_errors
#' @export
bootstrap <- function(data,
                      dep_col,
                      weight_vec,
                      sampling_method,
                      alpha,
                      jstar,
                      small,
                      trunc,
                      algorithm,
                      draw_weights,
                      var_names,
                      subsample_percent,
                      ...) {
  resample_qs(data, dep_col, weight_vec,
           sampling_method = "bootstrapRows", draw_weights = F,
           alpha, jstar, small, trunc, algorithm, subsample_percent = 1,
           var_names = var_names, ...)
}

#' @rdname standard_errors
subsample <- function(data,
                      dep_col,
                      sampling_method,
                      alpha,
                      jstar,
                      small,
                      trunc,
                      algorithm,
                      subsample_percent,
                      var_names,
                      draw_weights,
                      weight_vec,
                      ...) {
  resample_qs(data, dep_col, weight_vec,
           sampling_method = "subsampleRows", draw_weights = F,
           alpha, jstar, small, trunc, algorithm, subsample_percent = subsample_percent,
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
#' @param data Regression specification matrix.
#' @param dep_col Column of response variable.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param subsample_percent percentage of clusters to sample, must be greater than 0 but less than or equal to 1
#' @param draw_weights Boolean value; if true, draw a vector of exponential
#'                 weights to use in subsample
#' @param num_bs Number of subsample draws (must be greater than 1).
#' @param parallel whether to run in parallel or not
#' @param small Minimum size of residuals for computational accuracy.
#' @param trunc Boolean value; if true, replace those dependent values less than small with small itself;
#'         else, only use rows with residuals greater than small
#' @param cluster_matrix Matrix of cluster variables, as returned by a model formula
#' @param se_method Method to use for standard errors, either "weighted_bootstrap",
#' "subsample", "bootstrap" or "resample_qs" along with a specified subsampling method
#' (one of "leaveRows", "subsampleRows", or "bootstrapRows") and a subsampling percent.
#' @param sampling_method One of "leaveRows", "subsampleRows", or "bootstrapRows".
#' leaveRows doesn't resample rows at all. subsampleRows samples without replacement
#' given some percentage of the data (specified via subsample_percent), and boostrapRows
#' samples with replacement.
#' @param weight_vec vector of same length and order as dependent column, to be used as weights for estimation
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
standard_errors <- function(dep_col,
                            data,
                            cluster_matrix,
                            algorithm,
                            se_method,
                            var_names,
                            alpha,
                            jstar,
                            num_bs = 100,
                            parallel = F,
                            small = 1e-6,
                            trunc = FALSE,
                            weight_vec = NULL,
                            subsample_percent,
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
    cluster_index <- rep(1, nrow(data))
  } else {
    cluster_index <- get_strata(cluster_matrix)
  }

  num_clusters <- length(unique(cluster_index))


  fits <- future.apply::future_lapply(1:num_bs,
                                      future.seed = seed,
                            purrr::quietly(
                            function(x, dep_col,
                                     data,
                                     cluster_index,
                                     num_clusters,
                                     algorithm,
                                     se_method,
                                     var_names,
                                     alpha,
                                     jstar,
                                     small,
                                     trunc,
                                     weight_vec,
                                     subsample_percent,
                                     ...) {
        clustered_fits <- list()
        for(j in 1:num_clusters) {
          clustered_fits[[j]] <-
            do.call(se_method,
                  args = list(data = data[which(cluster_index == j),],
                              dep_col = dep_col[which(cluster_index == j)],
                              weight_vec = weight_vec[which(cluster_index == j)],
                              alpha = alpha, jstar = jstar, small = small,
                              trunc = trunc, algorithm = algorithm,
                              subsample_percent = subsample_percent,
                              var_names = var_names,
                              ...))
        }
        purrr::pmap(clustered_fits, rbind)
  }), data = data, se_method = se_method, cluster_index = cluster_index,
  dep_col = dep_col, weight_vec = weight_vec, alpha = alpha,
  jstar = jstar, var_names = var_names, small = small, trunc = trunc, algorithm = algorithm,
  subsample_percent = subsample_percent, num_clusters = num_clusters,  ...)

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
