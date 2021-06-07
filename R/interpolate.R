#' For sequencing along probabilities
#' @param from where to start
#' @param to where to end
#' @param by how far should each step be
seq_quant <- function(from, to, by) {
  s = seq(from, to, by)
  subset(s,s > 0 & s < 1)
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
    quantiles <- matrix(colMeans(object$quantreg_fit$quantiles), nrow = 1)
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
distributional_effects.matrix <- function(object, tails, alphas, ...) {

  params <- q_spline_R(quantiles = object,
                       alphas = alphas, tails = tails)

  pdf <- function(x) {
    if(length(x) > 1 & length(x) != nrow(object)) {
      stop("Vector length must either be 1 or equal to number of observations in new data.")
    }
    if(length(x) == 1) {
      x = rep(x, nrow(object))
    }
      splint_R(y = x,
               quantiles = object,
               alphas = alphas,
               y2 = params$y2, tail_param_l = params$tail_param_l,
               tail_param_u = params$tail_param_u,
               tails = tails, q_shift = params$q_shift,
               q_stretch = params$q_stretch,
               distn = "p")
  }

  cdf <- function(x) {
    if(length(x) > 1 & length(x) != nrow(object)) {
      stop("Vector length must either be 1 or equal to number of observations in new data.")
    }
    if(length(x) == 1) {
      x = rep(x, nrow(object))
    }
    splint_R(y = x,
             quantiles = object,
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

    if(length(x) > 1 & length(x) != nrow(object)) {
      stop("Vector length must either be 1 or equal to number of observations in new data.")
    }
    if(length(x) == 1) {
      x = rep(x, nrow(object))
    }
        splint_R(y = x,
                 quantiles = object,
                 alphas = alphas,
                 y2 = params$y2, tail_param_l = params$tail_param_l,
                 tail_param_u = params$tail_param_u,
                 tails = tails, q_shift = params$q_shift,
                 q_stretch = params$q_stretch,
                 distn = "q")
  }

  r <- function(n) {
    n_rows = nrow(object)
    out <- matrix(nrow = nrow(object), ncol = n)
    n_rows = nrow(object)

    for(i in 1:n) {
      x <- runif(n_rows, 0, 1)

      out[,i] = splint_R(y = x,
               quantiles = object,
               alphas = alphas,
               y2 = params$y2, tail_param_l = params$tail_param_l,
               tail_param_u = params$tail_param_u,
               tails = tails, q_shift = params$q_shift,
               q_stretch = params$q_stretch,
               distn = "q")
    }
    if(nrow(out) == 1) {
      out <- as.vector(out)
    }
    out
  }

  structure(list(
    pdf = pdf,
    cdf = cdf,
    q = q,
    r = r
  ), class = "distributional_effects")
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
    low <- min(x$q(tail_level))
    high <- max(x$q(1 - tail_level))
  }

  x_min <- fun(low)
  if(length(x_min) == 1) {
    vectorized_fun <- function(x) {
      sapply(x, fun)
    }
    x <- data.frame(x = c(low, high))
    ggplot2::ggplot(x, ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = vectorized_fun)
  } else {
    vectorized_fun <- function(x) {
      lapply(x, fun)
    }
    x = seq(low, high, by = (high - low) / 1000)
    y_mat = vectorized_fun(x)
    y_mat = do.call("cbind", y_mat)
    g <- ggplot2::ggplot(data.frame(x = x), ggplot2::aes(x = x))
    for(i in 1:nrow(y_mat)) {
      temp_data = data.frame(x = x, y = y_mat[i,])
      g <- g + ggplot2::geom_path(data = temp_data, ggplot2::aes(x = x, y = y), alpha = 0.8)
    }
    g
  }

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

  }

  interp <- future.apply::future_lapply(l, f, ...)

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
