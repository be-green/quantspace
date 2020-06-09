
#' Generate warning messages for regression model
#' @param model result of call to quantRegress
printWarnings = function(model){

  if(model$ierr != 0){
    warning(paste('Quantile Regresion ran into warning', model$ierr, ' and had ',model$it,' iterations.'))
    if(model$it < 12){
      warning('Iterations failed to go past threshold')
    }
  }
}
