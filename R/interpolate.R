#' For sequencing along probabilities
#' @param from where to start
#' @param to where to end
#' @param by how far should each step be
seq_quant <- function(from, to, by) {
  s = seq(from, to, by)
  subset(s,s > 0 & s < 1)
}

#' Simple vectorized version of splint_R
#' @param y passed to splint_R
#' @param ... other things passed to splint_R
#' @details This is an internal helper function
#' @importFrom future.apply future_sapply
vec_splint_R <- function(y, ...) {
  future.apply::future_sapply(y, splint_R, ...)
}

#' Interpolate a single set of quantiles
#' @param quantile matrix of fitted quantiles to interpolate
#' @param alphas values those quantiles are fit to
#' @param x values along which to evaluate the CDF, PDF, or quantile function
#' @param distn distribution to use, one of "q", "c", or "p"
#' @param tails what distribution to use for fitting tail distributions, either "gaussian"
#' or "exponential"
#' @details Used as a helper function for interpolate_quantiles
interpolate_one_row <- function(quantile, alphas, x, distn = "q", tails = "gaussian") {
  quantile <- matrix(quantile, nrow = 1)
  params <- q_spline_R(quantiles = quantile,
                       alphas = alphas, tails = tails)

  vec_splint_R(y = x,
               quantiles = quantile,
               alphas = alphas,
               y2 = params$y2, tail_param_l = params$tail_param_l,
               tail_param_u = params$tail_param_u,
               tails = tails, q_shift = params$q_shift,
               q_stretch = params$q_stretch,
               distn = distn)
}

#' Calculate distributional effects
#' @param object fit of class qs or matrix of fitted quantiles
#' @param newdata new data to predict distributional outcomes for
#' @param tails one of "gaussian" or "exponential"
#' @param alphas which quantiles these were fitted at
#' @param ... other parameters to pass
#' @details The arguments alphas and quantiles are automatically handled
#' if you pass an object of class "qs". If you don't pass new data to be
#' predicted on, it will assume that you want to calculate distributional
#' effects at the average of the data. This varies because of the non-linear
#' model for the quantile process
#' @export
distributional_effects <- function(object, tails = "gaussian", ...) {
  UseMethod("distributional_effects")
}

#' @rdname distributional_effects
#' @importFrom stats predict
#' @export
distributional_effects.qs <- function(object, tails = "gaussian",  newdata = NULL,  ...) {
  if(is.null(newdata)) {
    quantiles <- colMeans(object$quantreg_fit$quantiles)
  } else {
    quantiles <- predict(object, newdata = newdata)
  }

  alphas <- object$specs$alpha
  distributional_effects(object = quantiles, tails = tails, alphas = alphas)
}

#' helper function, borrowed from foreach
#' @param a value,
#' @param ... existing list
defcombine <- function(a, ...) {
  c(list(a), list(...))
}

#' helper function, collapse using correct method
#' @param x value,
#' @param l list
collapse_correctly <- function(x, l) {
  if(length(x) > 1) {
    do.call("rbind", l)
  } else {
    do.call("c", l)
  }
}

#' @rdname distributional_effects
#' @importFrom stats runif
#' @export
distributional_effects.numeric <- function(object, tails, alphas, ...) {

  quantiles <- matrix(object, nrow = 1)
  params <- q_spline_R(quantiles = quantiles,
                       alphas = alphas, tails = tails)

  pdf <- function(x) {
    vec_splint_R(y = x,
               quantiles = quantiles,
               alphas = alphas,
               y2 = params$y2, tail_param_l = params$tail_param_l,
               tail_param_u = params$tail_param_u,
               tails = tails, q_shift = params$q_shift,
               q_stretch = params$q_stretch,
               distn = "p")
  }

  cdf <- function(x) {
    vec_splint_R(y = x,
                 quantiles = quantiles,
                 alphas = alphas,
                 alphas = alphas,
                 y2 = params$y2, tail_param_l = params$tail_param_l,
                 tail_param_u = params$tail_param_u,
                 tails = tails, q_shift = params$q_shift,
                 q_stretch = params$q_stretch,
                 distn = "c")
  }

  q <- function(x) {
    if(x > 1 | x < 0) {
      stop("Probability must be between 0 & 1")
    }
    vec_splint_R(y = x,
                 quantiles = quantiles,
                 alphas = alphas,
                 y2 = params$y2, tail_param_l = params$tail_param_l,
                 tail_param_u = params$tail_param_u,
                 tails = tails, q_shift = params$q_shift,
                 q_stretch = params$q_stretch,
                 distn = "q")
  }

  r <- function(n) {

    x <- runif(n, 0, 1)

    vec_splint_R(y = x,
                 quantiles = quantiles,
                 alphas = alphas,
                 y2 = params$y2, tail_param_l = params$tail_param_l,
                 tail_param_u = params$tail_param_u,
                 tails = tails, q_shift = params$q_shift,
                 q_stretch = params$q_stretch,
                 distn = "q")
  }



  structure(list(
    pdf = pdf,
    cdf = cdf,
    q = q,
    r = r
  ), class = "distributional_effects")
}


#' @rdname distributional_effects
#' @export
distributional_effects.matrix <- function(object, tails = "gaussian", alphas, ...) {

  fitted <- map_rows_parallel(object, f = distributional_effects,
                              alphas = alphas, tails = tails, ..., collapse = "list")

  pdf <- function(x) {
    l <- lapply(fitted, function(f) f$pdf(x))
    if(length(x) > 1) {
      do.call("rbind", l)
    } else {
      do.call("c", l)
    }
  }

  cdf <- function(x) {
    l <- lapply(fitted, function(f) f$cdf(x))
    if(length(x) > 1) {
      do.call("rbind", l)
    } else {
      do.call("c", l)
    }
  }

  q <- function(x) {
    l <- lapply(fitted, function(f) f$q(x))
    if(length(x) > 1) {
      do.call("rbind", l)
    } else {
      do.call("c", l)
    }
  }

  r <- function(n) {
    l <- lapply(fitted, function(f) f$r(n))

    if(n > 1) {
      do.call("rbind", l)
    } else {
      do.call("c", l)
    }
  }


  structure(list(
    fitted = fitted,
    pdf = pdf,
    cdf = cdf,
    q = q,
    r = r
  ), class = "distributional_effects_list")
}

#' Visualize distributional effects
#' @param x estimate of distributional effects
#' @param what what to plot, one of "pdf", "cdf", "quantiles"
#' @param tail_level probability level to stop plotting tail
#' @param ... other arguments, ignored for now
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 stat_function
#' @export
plot.distributional_effects <- function(x,
                                        what = "pdf",
                                        tail_level = 0.01,
                                        ...) {

  assertthat::assert_that(0 < tail_level)
  assertthat::assert_that(1 > tail_level)
  assertthat::assert_that(length(tail_level) == 1)
  assertthat::assert_that(is.numeric(tail_level))
  assertthat::assert_that(is.character(what))

  if(what == "quantiles") {
    what = "q"
  }

  fun <- x[[what]]

  if(what == "q") {
    low = tail_level
    high = 1 - tail_level
  } else {
    low <- x$q(tail_level)
    high <- x$q(1 - tail_level)
  }

  x <- data.frame(x = c(low, high))
  ggplot2::ggplot(x, ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = fun)
}

#' Visualize distributional effects
#' @param x estimate of distributional effects
#' @param what what to plot, one of "pdf", "cdf", "quantiles"
#' @param tail_level probability level to stop plotting tail
#' @param ... other arguments, ignored for now
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 stat_function
#' @export
plot.distributional_effects_list <- function(x,
                                        what = "cdf",
                                        tail_level = 0.01,
                                        ...) {

  assertthat::assert_that(0 < tail_level)
  assertthat::assert_that(1 > tail_level)
  assertthat::assert_that(length(tail_level) == 1)
  assertthat::assert_that(is.numeric(tail_level))
  assertthat::assert_that(is.character(what))


  if(what == "quantiles") {
    what = "q"
  }

  g <- ggplot2::ggplot(data.frame())

  distributional_effects_list <- x$fitted

  N = length(distributional_effects_list)

  for(i in 1:N) {

    de <- distributional_effects_list[[i]]
    fun <- de[[what]]

    if(what == "q") {
      low = tail_level
      high = 1 - tail_level
    } else {
      low <- de$q(tail_level)
      high <- de$q(1 - tail_level)
    }

    x <- data.frame(x = c(low, high))
    g <- g + ggplot2::stat_function(data = x,
                                    ggplot2::aes(x = x), fun = fun, alpha = 0.8)
  }

  g
}

#' Map a function along a list in parallel
#' @param l list whose items are passed as data to f
#' @param f function to map along rows
#' @param ... parameters passed to function
#' @param parallel T/F, whether to operate in parallel or in sequence
#' @param ncores number of cores to use
#' @param thresh required length of list to use parallel processing
#' @param collapse what function to use when collapsing arguments
#' @details only works if the function's dependencies are completely
#' contained in "quantspace" package
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future.apply future_lapply
map_parallel <- function(l, f, ..., parallel = T,
                         ncores = getCores(), thresh = 20,
                         collapse = "list") {
  if(parallel & (length(l) > thresh)) {

    ncores <- getCores()
    setCores(ncores)

    interp <- future.apply::future_lapply(l, f, ...)

  } else {

    old_plan <- future::plan()
    future::plan(future::sequential)

    interp <- list()
    for(i in 1:length(l)) {
      interp[[i]] <- f(l[[i]], ...)
    }

    future::plan(old_plan)

  }
  do.call(collapse,interp)
}

#' Map a function along rows of a matrix or data.frame
#' @param mat matrix or data.frame whose rows are passed as data to f
#' @param f function to map along rows
#' @param ... parameters passed to function
#' @param parallel T/F, whether to operate in parallel or in sequence
#' @param collapse function to use when collapsing list of objects
#' @param ncores number of cores to use
#' @param row_thresh required number of rows to use parallel processing
#' @importFrom future.apply future_apply
#' @importFrom future nbrOfWorkers
#' @importFrom future sequential
#' @importFrom future plan
#' @details only works if the function's dependencies are completely
#' contained in "quantspace" package
map_rows_parallel <- function(mat, f, ..., parallel = T,
                              ncores = getCores(), row_thresh = 20,
                              collapse = "rbind") {
  if(parallel & (nrow(mat) > row_thresh)) {

    if(ncores != future::nbrOfWorkers()) {
      makePlan(ncores)
    }

    interp <- future.apply::future_apply(mat, MARGIN = 1, FUN = f,
                                         simplify = F, ...)

  } else {

    old_plan <- future::plan()
    future::plan(future::sequential)

    interp <- list()
    for(i in 1:nrow(mat)) {
      interp[[i]] <- f(mat[i,], ...)
    }
    interp <- (do.call(collapse,interp))
    future::plan(old_plan)
  }

  do.call(collapse,interp)

}

#' Interpolate quantiles and return a cumulative distribution function
#' @param object either matrix of quantiles or a fitted model with a relevant
#' method to dispatch
#' @param ... other arguments passed to interpolate quantiles
#' @param alphas values those quantiles are fit to
#' @param grid grid along which to evaluate the CDF
#' @param parallel whether to work in parallel
#' @param ncores number of cores to use
#' @param distn what to return, q for quantile, c for cdf, p for pdf
#' @export
interpolate_quantiles <- function(object, ...) {
  UseMethod("interpolate_quantiles")
}

#' @rdname interpolate_quantiles
#' @param newdata New data to use when calculating quantiles to interpolate
#' @param grid grid along which to evaluate the CDF
#' @param parallel whether to work in parallel
#' @param ncores number of cores to use
#' @param tails what distribution to use when fitting the tails, either "gaussian" or "exponential"
#' @param distn what to return, q for quantile, c for cdf, p for pdf
#' @export
interpolate_quantiles.qs <- function(object, newdata = NULL,
                                     grid = seq_quant(0, 1, by = 0.01),
                                     parallel = T, ncores = getCores(),
                                     row_thresh = 20, tails = "gaussian",
                                     distn = "q", ...) {
  quantiles <- predict(object, newdata = newdata)
  alphas <- object$specs$alpha

  interpolate_quantiles(quantiles, alphas, grid, parallel, ncores,
                        row_thresh = 20, tails = tails,
                        distn = distn)
}

#' @rdname interpolate_quantiles
#' @param alphas values those quantiles are fit to
#' @param grid grid along which to evaluate the CDF
#' @param parallel whether to work in parallel
#' @param ncores number of cores to use
#' @param distn what to return, q for quantile, c for cdf, p for pdf
#' @param row_thresh required minimum number of observations to use parallel processing
#' @export
interpolate_quantiles.matrix <- function(object, alphas,
                                         grid = seq_quant(0, 1, by = 0.01),
                                  parallel = T, ncores = getCores(),
                                  row_thresh = 100, tails = "gaussian",
                                  distn = "q", ...) {

  interp <- map_rows_parallel(
    object, f = interpolate_one_row,
                             alphas = alphas,
                             x = grid,
                             distn = "q",
    parallel = parallel, ncores = getCores(),
    row_thresh = row_thresh,
    tails = tails,
    distn = distn
  )

  colnames(interp) <- grid
  rownames(interp) <- NULL
  return(interp)
}

