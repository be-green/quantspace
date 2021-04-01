#' Set bootstrap cores on load + welcome messages
#' @param libname legacy requirement
#' @param pkgname legacy requirement
#' @importFrom parallel detectCores
.onAttach <- function(libname, pkgname) {

  ncores = getOption('qs.cores')
  if(is.null(ncores)) {
    ncores = round(parallel::detectCores()/2)
  }

  setCores(ncores)

  packageStartupMessage("Loaded quantspace v0.1, using ", ncores,
                        " cores for bootstrap sampling (see ?getCores).\n",
                        "Bug reports: github.com/be-green/quantspace/issues")
}
