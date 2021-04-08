#' Set bootstrap cores on load + welcome messages
#' @param libname legacy requirement
#' @param pkgname legacy requirement
#' @importFrom parallel detectCores
#' @importFrom future sequential
#' @importFrom future plan
#' @importFrom future availableCores
.onAttach <- function(libname, pkgname) {

  currently_sequential = is_sequential()

  ncores = getOption('qs.cores')
  if(is.null(ncores) & currently_sequential) {
    ncores = future::availableCores()
  }

  if(ncores >= 1) {
    setCores(ncores)
  }

  packageStartupMessage("Loaded quantspace v0.1, using ", ncores,
                        " cores for bootstrap sampling (see ?getCores).\n",
                        "Bug reports: github.com/be-green/quantspace/issues")
}

is_sequential <- function() {
  "sequential" %in% setdiff(class(future::plan()), c("FutureStrategy", "tweaked",
                            "function"))
}
