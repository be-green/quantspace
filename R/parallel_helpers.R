
#' Get user defined cores
#' @rdname getCores
#' @export
getCores <- function() {
  mc_cores <- getOption("mc.cores")
  if(is.null(mc_cores)) {
    1
  } else {
    max(c(1, mc_cores))
  }
}

#' @rdname getCores
#' @param ncores number of cores to use for bootstrapping
#' @importFrom assertthat assert_that
#' @export
#' @details `getCores()` and `setCores()` determine the number of cores
#' used in the boostrap procedure used to generate standard errors for the
#' quantile spacings estimator. By default, this is set to half of available
#' cores, or whatever is found in `getOption('mc.cores')`. This package uses
#' the `future` backend for parallelization, so if you would like to specify
#' a custom plan for futures use `future::plan()` after loading the package.
#' A better interface for custom plans is on the roadmap for the package but
#' isn't currently implemented.
setCores <- function(ncores) {
  assertthat::assert_that(is.numeric(ncores))
  assertthat::assert_that(ncores > 0)
  options(mc.cores = ncores)
  makePlan(ncores)
}

#' Make a plan for the future parallel backend
#' @param ncores number of cores
#' @details This function gets the number of cores set by user
#' and subsequently makes a plan that is operating system dependent.
#' If it is on windows, it uses multisession since forks are not an option.
#' Otherwise it uses multicore. If you would like to pass a custom plan.
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @importFrom future nbrOfWorkers
makePlan <- function(ncores) {
  assertthat::assert_that(is.numeric(ncores))
  assertthat::assert_that(ncores > 0)
  if(ncores == 1) {
    future::plan(future::sequential)
  } else if(ncores == future::nbrOfWorkers()) {
    # don't change anything
  } else if(future::supportsMulticore()) {
    future::plan(future::multicore, workers = ncores)
  } else {
    future::plan(future::multisession,workers = ncores)
  }
}
