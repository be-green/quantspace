
#' Draw a subsample of clusters and optionally assigns random exponential weights
#' @param cluster_col Column vector with list of cluster indices
#' @return N_cluster x 2 vector with starting
#'                   and ending positions of each cluster.
#'                   NOTE that this requires that data are
#'                   sorted in this dimension
getClusterIndices = function(cluster_col) {

  # TODO: probably should add some logic to check sorting...

  uniqClusters<-unique(as.vector(cluster_col))
  cluster_indices<-matrix(0,length(uniqClusters),2)
  for (jj in 1:length(uniqClusters)){
    index_find<-which(cluster_col==uniqClusters[jj])
    cluster_indices[jj, 1] = index_find[1]
    cluster_indices[jj, 2] = index_find[length(index_find)]
  }

  return(cluster_indices)

}

#' Draws a subsample of clusters and (optionally) assigns
#' random, exponential weights to each cluster
#' @param cluster_indices N_cluster x 2 vector with starting
#' and ending positions of each cluster.
#' NOTE that this requires that data are
#' sorted in this dimension
#' @param stratum_indices N_cluster x 1 vector with stratum indices
#' (factor variables). This enables the
#' researcher to stratify (e.g., on cluster size).
#' Set to NULL to take a random subsample.
#' @param M percentage of clusters to sample.  Must be between 0 and 1 (inclusive).
#' @param draw_weights If true, draw a vector of exponential weights to use in subsample
#' @return list of
#' subsample_rows: Boolean vector of selected rows of sampled clusters
#' subsample_weights: either NULL, or vector of cluster weights
#' @importFrom stats rexp
clusterSample = function(cluster_indices,
                         stratum_indices = NULL,
                         M = 1,
                         draw_weights = FALSE) {

  # TODO: probably should add some logic to validate inputs here...

  # checking conditions on M
  if(M > 1 | 0 >= M){
    stop("Subsampling percentage must be between 0 and 1")
  }

  # draw an M% stratified sample of indices
  if (!is.null(stratum_indices)){
    sampled_indices <- do.call(rbind,
                               lapply(split(as.data.frame(cluster_indices), stratum_indices),
                                      function(x) x[sample(nrow(x), floor(nrow(x)*M)),]))
  } else {
    sampled_indices <- as.data.frame(cluster_indices[
      sample(nrow(cluster_indices),
             floor(nrow(cluster_indices)*M) ),
    ])
  }

  # now, re-sort the vector
  sampled_indices <- sampled_indices[order(sampled_indices[,1]),]
  n_sampled_rows <- sum((sampled_indices$V2-sampled_indices[,1]+1))
  # types seem to get messed up easily, so I'll just initialize with the sample function
  subsample_rows <- matrix(FALSE,max(cluster_indices[,2]),1)

  # drawing weights (if applicable)
  if (!draw_weights){
    subsample_weights <- NULL
  } else {
    cluster_weights <- pmax(stats::rexp(dim(sampled_indices)[1]),5e-3)
    subsample_weights <- matrix(0,n_sampled_rows,1)
  }

  # looping over sampled clusters to populate output vectors
  row_s_ct <- 0
  for (j in 1:dim(sampled_indices)[1] ){
    sidx <- sampled_indices[j,1]
    eidx <- sampled_indices[j,2]

    subsample_rows[sidx:eidx] <- TRUE

    if (draw_weights){
      subsample_weights[(row_s_ct+1):(row_s_ct+(eidx-sidx+1)),] <- cluster_weights[j]
    }
    row_s_ct <- row_s_ct + (eidx-sidx+1)
  }
  subsample_rows <- which(subsample_rows)

  return(list('subsample_rows' = subsample_rows,
              'subsample_weights' = subsample_weights))

}

#' Local do parallel function
#' @inheritParams foreach::`%dopar%`
#' @importFrom foreach `%dopar%`
`%dopar%` <- foreach::`%dopar%`

#' Get user defined cores
#' @rdname getCores
#' @export
getCores <- function() {
  mc_cores <- getOption("qs.cores")
  if(is.null(mc_cores)) {
    1
  } else {
    max(c(1, mc_cores))
  }
}

#' @rdname getCores
#' @param ncores number of cores to use for bootstrapping
#' @export
#' @details `getCores()` and `setCores()` determine the number of cores
#' used in the boostrap procedure used to generate standard errors for the
#' quantile spacings estimator. By default, this is set to half of available
#' cores, or whatever is found in `getOption('qs.cores')`.
setCores <- function(ncores) {
  assertthat::assert_that(is.numeric(ncores))
  assertthat::assert_that(ncores > 0)
  options(qs.cores = ncores)
}

#' Computes standard errors for the quantile regression spacing method using
#' subsampling.
#' @param data Regression specification matrix.
#' @param dep_col Column of response variable.
#' @param var_names RHS regression variable names.
#' @param alpha Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param cluster_indices N_cluster x 2 vector with starting
#'                   and ending positions of each cluster.
#'                   NOTE that this requires that data are
#'                   sorted in this dimension
#' @param stratum_indices N_cluster x 1 vector with stratum indices
#'                   (factor variables). This enables the
#'                   researcher to stratify (e.g., on cluster size)
#'                   Set to NULL to take a random subsample
#' @param M percentage of clusters to sample, must be greater than 0 but less than or equal to 1
#' @param draw_weights Boolean value; if true, draw a vector of exponential
#'                 weights to use in subsample
#' @param num_bs Number of subsample draws (must be greater than 1).
#' @param parallel whether to run in parallel or not
#' @param num_cores number of cores to use, defaults to option set by `options(mc.cores)` if not specified
#' @param small Minimum size of residuals for computational accuracy.
#' @param trunc Boolean value; if true, replace those dependent values less than small with small itself;
#'         else, only use rows with residuals greater than small
#' @param start_model Starting values for regression's optimization.
#' @param weight_vec vector of same length and order as dependent column, to be used as weights for estimation
#'              (note, if draw weights is set to TRUE, this variable will be the element-wise product
#'              of itself and a random vector of weights)
#' @param algorithm function which is actually used to fit each quantile regression
#' @param ... other arguments passed to quantile fitting function
#' @return
#'   list of cov: num_betas x num_betas covariance matrix using
#'                 bootstrapped subsampling covariances
#'           ols_cov: VCV matrix from OLS estimates
#'           warnings: num_bs x p matrix of warning messages
#'                      produced in each bootstrap's quantile regression call
#'          iter: num_bs x p matrix of iterations ran by each
#'                 bootstrapped quantile regression call
#'          OLS: OLS point estimates
#'          counts: number of observations in each subsample
#'          coef_boot: full set of boostrapped quantile spacing coefficients
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParWorkers
#' @importFrom foreach getDoParVersion
#' @importFrom foreach getDoParName
#' @importFrom foreach foreach
#' @importFrom foreach `%dopar%`
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom methods is
#' @importFrom stats rexp
#' @importFrom stats cov
#' @importFrom purrr pmap
#' @importFrom stats dnorm
#' @export
subsampleStandardErrors = function(
  dep_col,
  data,
  algorithm,
  var_names,
  alpha,
  jstar,
  cluster_indices = NULL,
  stratum_indices = NULL,
  M = 1,
  draw_weights = TRUE,
  num_bs = 100, parallel = F, num_cores = getCores(), small = 1e-6,
  trunc = FALSE, start_model, weight_vec = NULL,
  ...) {

  # check to see if regression matrix is sparse. If not, then turn into CSR matrix
  if(!methods::is(data, 'matrix.csr')){
    data = denseMatrixToSparse(data)
  }

  # checking conditions on M
  if(M > 1 | 0 >= M){
    stop("Subsampling percentage must be between 0 and 1")
  }


  if(parallel){

    cl = parallel::makeCluster(num_cores, outfile = '')

    doParallel::registerDoParallel(cl)

    fit = foreach::foreach(
      bs = 1:num_bs,
      .packages = c('quantreg'),
      .export = c('quantRegSpacing', 'rq.fit.sfn_start_val','findRedundantCols',
                  'ensureSpecFullRank', 'getRank', 'addMissingSpecColumns', 'clusterSample',
                  'rho', 'getColNums', 'regressResiduals', 'printWarnings'),
      .errorhandling = 'stop' #, .options.snow = opts
    ) %dopar% {
      # if no clustering is specified, draw sample here.
      # Need to provide a cluster_id variable to stratify
      if (is.null(cluster_indices)) {
        rows = sample(dim(data)[1],floor(M*dim(data)[1]))
        if (draw_weights) {
          if(any(is.null(weight_vec))){
            rand_weight_vec <- pmax(stats::rexp(length(rows)),5e-3)
          } else {
            rand_weight_vec <- pmax(rexp(length(rows)),5e-3) * weight_vec[rows]
          }
        } else {
          rand_weight_vec <- weight_vec[rows]
        }
      } else {
        # call clusterSample here to generate clustered standard errors
        subsample_outputs <- clusterSample(cluster_indices = cluster_indices,
                                           stratum_indices = stratum_indices,
                                           M = M,
                                           draw_weights = draw_weights)
        rows <- subsample_outputs$subsample_rows
        if(any(is.null(weight_vec))) {
          rand_weight_vec <- subsample_outputs$subsample_weights
        } else{
          if(draw_weights) {
            rand_weight_vec <- subsample_outputs$subsample_weights * weight_vec[rows]
          }else {
            rand_weight_vec <- weight_vec[rows]
          }
        }
      }
      if (!exists("rand_weight_vec")) {
        rand_weight_vec <- NULL
      }

      cur_fit <- quantRegSpacing(
        data = data[rows,],
        dep_col = dep_col[rows],
        var_names = var_names,
        alpha = alpha,
        jstar = jstar,
        small = small,
        trunc = trunc,
        start_list = NA,
        weight_vec = rand_weight_vec,
        algorithm = algorithm,
        ...)

      return(cur_fit)
    }
    parallel::stopCluster(cl)
  } else {
    fit <- list()
    for(bs in 1:num_bs) {
      # if no clustering is specified, draw sample here.
      # Need to provide a cluster_id variable to stratify
      if (is.null(cluster_indices)) {
        rows = sample(dim(data)[1],floor(M*dim(data)[1]))
        if (draw_weights) {
          if(any(is.null(weight_vec))){
            rand_weight_vec <- pmax(stats::rexp(length(rows)),5e-3)
          } else {
            rand_weight_vec <- pmax(rexp(length(rows)),5e-3) * weight_vec[rows]
          }
        } else {
          rand_weight_vec <- weight_vec[rows]
        }
      } else {
        # call clusterSample here to generate clustered standard errors
        subsample_outputs <- clusterSample(cluster_indices = cluster_indices,
                                           stratum_indices = stratum_indices,
                                           M = M,
                                           draw_weights = draw_weights)
        rows <- subsample_outputs$subsample_rows
        if(any(is.null(weight_vec))) {
          rand_weight_vec <- subsample_outputs$subsample_weights
        } else{
          if(draw_weights) {
            rand_weight_vec <- subsample_outputs$subsample_weights * weight_vec[rows]
          } else {
            rand_weight_vec <- weight_vec[rows]
          }
        }
      }

      if (!exists("rand_weight_vec")) {
        rand_weight_vec <- NULL
      }

      cur_fit <- quantRegSpacing(
        data = data[rows,],
        dep_col = dep_col[rows],
        var_names = var_names,
        alpha = alpha,
        jstar = jstar,
        small = small,
        trunc = trunc,
        start_list = NA,
        weight_vec = rand_weight_vec,
        algorithm = algorithm,
        ...)

      fit[[bs]] <- cur_fit
    }
  }

  fit <- purrr::pmap(fit, rbind)
  quant_cov_mat <- stats::cov(fit$coef, use = "pairwise.complete.obs")
  quant_cov_mat <- quant_cov_mat * (M)

  return(list('quant_cov' = quant_cov_mat,
              'warnings' = fit$warnings,
              'iter' = fit$iter,
              'counts' = fit$counts,
              'coef_boot' = fit$coef,
              'M' = M))
}
