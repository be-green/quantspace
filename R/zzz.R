#' Set bootstrap cores on load + welcome messages
#' @param libname legacy requirement
#' @param pkgname legacy requirement
#' @importFrom parallel detectCores
#' @importFrom future sequential
#' @importFrom future plan
#' @importFrom future availableCores
.onAttach <- function(libname, pkgname) {

  ncores = getOption('mc.cores')
  if(is.null(ncores)) {
    ncores = pmax(round(future::availableCores()/2), 1)
    options(mc.cores = ncores)
  }

  packageStartupMessage("Loaded quantspace v0.1, using ", ncores,
                        " cores for bootstrap sampling (see ?getCores).\n",
                        "Bug reports: github.com/be-green/quantspace/issues")
}

