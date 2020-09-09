#' Get combinations of two factors
#' @param x first factor
#' @param y second factor
get_combos <- function(x, y) {
  suppressWarnings({levels(interaction(x, y))})
}

#' Calculate r^{k+1} given parameters
#' Follows Yu and Lin (2017)
#' @param rho relaxation parameter
#' @param u Lagrange multiplier term
calc_r <- function(rho, u, Y, Gamma, tau, One) {
  pmax((1 / rho) * u + Y - Gamma -(tau / rho) * One, 0 * One) - pmax(-(1 / rho) * u - Y + Gamma + ((tau - 1) / rho) * One, 0 * One)
}


#' Update gamma given an initial matrix
#' @param Gamma Parameter matrix
#' @param tune Tuning parameter to pass along
#' @param W Number of groups at a given time
#' @param nu Penalization parameter used in optimization
#' @param E Proportion of observed values at a given time (correction for imbalanced panel)
#' @param err_g_thresh error threshold for gamma calculation
#' @importFrom data.table setnames
#' @importFrom data.table as.data.table
update_gamma <- function(Gamma, tune, W, E, nu, err_g_thresh = 0.01) {
  err_g <- 1000
  while (err_g > err_g_thresh) {
    # Iterate as in the note
    S <- svd(as.matrix(Gamma - tune * W * (W * Gamma - E * W)))
    D <- as.vector(pmax(S$d - matrix(1, length(S$d), 1) * tune * nu / 2, 0))
    Gamma_New <-  S$u %*% diag(D, nrow(S$v)) %*% t(S$v)

    err_g <- (norm(as.matrix(Gamma)-Gamma_New)/(norm(as.matrix(Gamma))+1))
    Gamma <- data.table::as.data.table(Gamma_New)
    if(is.na(err_g)) {
      browser()
    }
    # print(Error_G)
  }
  data.table::setnames(Gamma, colnames(Gamma), as.character(1:ncol(Gamma)))
  Gamma
}

#' Replace all NA values in a data.table
#' @param DT data.table to put 0 values in (instead of NA values)
#' @details edits values by reference, so no need to return
#' or assign things. I know, not very functional programming,
#' but it is really fast.
#' credit: https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
replace_nas <- function(DT) {
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    data.table::set(DT,which(is.na(DT[[j]])),j,0)
}

#' Calculate quantile regression via admm estimation
#' @param Y outcome variable
#' @param X variables to estimate
#' @param vec vector of parameters to update
#' @param rho stochastic relaxation parameter
#' @param tau quantile to estimate
#' @param k penalization parameter
#' @param iterfun function for updating parameters given new errors
#' @param err_u_thresh error threshold for convergence
#' @param verbose whether to be chatty
#' @param ... other parameters to be passed to f
admm <- function(Y, X, est, vec, rho, tau, iterfun, k, err_u_thresh, verbose, ...) {

  err_u <- 1000
  iter <- 1

  DT <- data.table(Y, X, vec)

  # Big matrix of ones if needed
  One <- matrix(1, nrow = length(Y), ncol = 1)

  while(err_u > err_u_thresh) {

    if (iter == 1) {
      # Initial u
      DT[,u := rep(0, length(Y))]
    } else {

      err_u <- norm(as.matrix(rho * (DT$Y - DT$r - DT$vec)))/(sum(abs(DT$u)) + 1)
      DT[,u := u + rho * (Y - r - vec)]

      if(verbose) {
        message(paste0("Iteration: ", iter, "\n",
                       "Error: ", signif(err_u, 4)))
      }

    }

    # Append to data frame
    DT[,r := calc_r(rho, u, Y, vec, tau, One)]
    DT[,e := Y - r + u / rho]

    # update estimate and predictions
    est <- iterfun(e = DT$e, vec = DT$vec, est, ...)
    vec <- predict(est, DT, ...)
    DT[,vec := vec]
    iter = iter + 1

  }
  est
}


#' Iteration for admm algorithm for two-way fixed effects
#' @param DT data.table with all the model info
#' @param est previous estimate of Gamma
#' @param Year vector of year assignments
#' @param Group vector of group assignments
#' @param tune tuning parameter for updates
#' @param nu Penalization parameter used in optimization
#' @param W matrix or data.table of sqrt(N) by group
#' @param J_levels factor representing levels of grouping variable, unused here
#' @param T_levels factor representing levels of grouping variable, unused here
#' @param err_g_thresh threshold for convergence of gamma updates
gamma_admm_iter <- function(e, vec, est, W, Year, Group, tune, nu, J_levels, T_levels,
                            err_g_thresh) {

  if(anyNA(est) | anyNA(vec)) {
    browser()
  }

  e_dt <- data.table::data.table(e, Year, Group)

  Gamma <- data.table::as.data.table(est)

  # Ugly but fast
  # Would be nice to avoid this if possible
  E <- data.table::dcast(
    e_dt[,.(e1 = sum(e)/length(e)),
       by = c("Year", "Group")],
    ... ~ Year, value.var = "e1"
  )[
    ,Group := NULL
  ]
  replace_nas(E)

  # Update Gamma
  Gamma <- update_gamma(Gamma = Gamma,
                        tune = tune,
                        nu = nu,
                        W = W,
                        E = E,
                        err_g_thresh = err_g_thresh)

  structure(Gamma, class = c("two_way_fe_estimate", "data.table"))
}

#' Prediction function for two-way fixed effects given DT
#' @param est previous estmiate of two way quantile fixed effects
#' @param DT data.table that contains design matrix columns to allow prediction
#' @details Right now this is super fragile, expecting specifically
#' one column named "Group" and one named "Year". This is not meant to be viewed
#' externally to the package.
predict.two_way_fe_estimate <- function(est, DT, ...) {
  merge_gamma <- melt_gamma(est, J_levels = DT$Group, T_levels = DT$Year)
  DT <- merge_gamma[DT, on = c("Group", "Year")]
  DT$Gamma
}

#' Calculate two-way quantile fixed effects via admm
#' @param Y vector of outcomes
#' @param T_vec vector of time assignment
#' @param J_vec vector of group assignments
#' @param N number of assets
#' @param Time integer representing number of periods
#' @param J integer representing number of groups
#' @param tau quantile to calculate
#' @param rho stochastic relaxation parameter
#' @param k penalization parameter used in iteration, defaults to no penalization
#' @param err_u_thresh error threshold for u
#' @param err_g_thresh error threshold for g
#' @importFrom data.table setkey
#' @importFrom data.table melt
#' @importFrom data.table dcast
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table `:=`
#' @importFrom data.table setorder
#' @importFrom data.table setnames
#' @importFrom data.table copy
#' @import data.table
#' @export
calc_gamma <- function(Y,
                       T_vec,
                       J_vec,
                       N,
                       J,
                       Time,
                       Gamma = NULL,
                       tau = 0.5,
                       rho = 0.1,
                       k = 1,
                       err_u_thresh = 0.01,
                       err_g_thresh = 0.01,
                       verbose = F) {


  J_levels <- factor(unique(J_vec))
  T_levels <- factor(unique(T_vec))

  DT <- data.table::data.table(Y = Y,
                               Group = factor(J_vec, levels = J_levels),
                               Year = factor(T_vec, levels = T_levels))

  data.table::setorder(DT, Year, Group)

  NT <- length(Y)

  interactions <- data.table(Interactions = get_combos(T_levels,J_levels))

  DT[,Interactions := paste0(Year, ".", Group)]

  data.table::setkey(DT, "Interactions")
  data.table::setkey(interactions, "Interactions")


  nu <- k*max(sqrt(N), sqrt(Time))
  if(is.na(nu)) {
    browser()
  }

  # W is simply the counts of firms in each group
  W <- DT[,.(w = .N),
          by = c("Year", "Group")]

  # Set tuning parameter
  tune <- (1 / max(W$w))

  # One column per group
  W <- dcast(W, ... ~ Year, value.var = "w")[
    ,Group := NULL
  ]

  # replace by reference, no need for assignment
  replace_nas(W)

  # Need square roots, actually (see the notes)
  W <- sqrt(W)

  if(is.null(Gamma)) {
    vec <- rnorm(length(Y), mean = 0, sd = 1)
    Gamma <- data.table::as.data.table(matrix(rnorm(J * Time), J, Time))
  } else {
    vec <- predict(Gamma, DT[,.(Y, Group, Year)], J_levels, T_levels)
  }

  Gamma <- admm(Y, X = DT[,.(Group, Year)], est = Gamma, vec = vec, rho = rho,
                tau = tau, iterfun = gamma_admm_iter,
                k = k,
                err_u_thresh = err_u_thresh,
                verbose = verbose,
                W = W, Year = DT$Year, Group = DT$Group,
                J_levels = J_levels, T_levels = T_levels,
                tune = tune,
                err_g_thresh = err_g_thresh,
                nu = nu)

  Gamma
}

#' Reformat Gamma matrix of quantile estimates for fixed effects
#' @param Gamma matrix of coefficients for fixed effects
#' @param J_levels levels of grouping factor
#' @param T_levels levels of time factor
#' @importFrom data.table melt
#' @importFrom data.table copy
#' @importFrom data.table `:=`
melt_gamma <- function(Gamma, J_levels, T_levels) {
  data.table::melt(data.table::copy(data.table::as.data.table(Gamma))[,Group := factor(1:.N, labels = levels(J_levels))],
                   id.vars = "Group", variable.factor = F, variable.name = "Year",
                   value.name = "Gamma")[,
                                         Year := factor(as.integer(gsub("[^0-9]", "",
                                                                        Year)),
                                                        labels = levels(T_levels))]
}

#' Reformat Gamma vector of quantile estimates for fixed effects
#' @param Gamma vector of coefficients for fixed effects
#' @param J_levels levels of grouping factor
#' @param T_levels levels of time factor
#' @importFrom data.table melt
#' @importFrom data.table copy
#' @importFrom data.table `:=`
cast_gamma <- function(Gamma, J_vec, T_vec) {
  data.table::dcast(unique(data.table::data.table(Gamma,Group = J_vec, Year = T_vec)),
                    ...~Year, value.var = "Gamma")[
                      ,Group := NULL
                    ][]
}

#' Iterate over K until getting a low enough rank for Gamma
#' @param Y vector of outcomes
#' @param T_vec vector of time assignment
#' @param J_vec vector of group assignments
#' @param N number of assets
#' @param Time integer representing number of periods
#' @param J integer representing number of groups
#' @param tau quantile to calculate
#' @param rho stochastic relaxation parameter
#' @param k_init initial penalization parameter used in iteration, defaults to 1
#' @param err_u_thresh error threshold for u
#' @param err_g_thresh error threshold for g
#' @importFrom data.table setkey
#' @importFrom data.table melt
#' @importFrom data.table dcast
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table `:=`
#' @importFrom data.table setorder
#' @importFrom data.table setnames
#' @importFrom data.table copy
#' @import data.table
#' @importFrom Matrix rankMatrix
#' @export
calc_reduce_rank_gamma <- function(Y,
                              T_vec,
                              J_vec,
                              N,
                              J,
                              Time,
                              Gamma = NULL,
                              tau = 0.5,
                              rho = 0.1,
                              init_k_high = 30,
                              init_k_low = 1,
                              desired_dim = NULL,
                              err_u_thresh = 0.01,
                              err_g_thresh = 0.01,
                              verbose = F,
                              maxiter = 100,
                              parallel = T) {

  if(is.null(desired_dim)) {
    stop("Must specify target dimensionality for Gamma.")
  }



  # initial window for k
  lk <- init_k_low
  hk <- init_k_high

  high_gamma <- NULL
  low_gamma <- NULL

  stop <- 0
  RANK <- "Unknown"

  if(parallel == T) {
    future::plan(future::multiprocess(workers = 2))
  }

  lr <- vector(mode = "numeric")
  hr <- vector(mode = "numeric")


  i = 1

  cur_hr <- NA
  cur_lr <- NA

  while(is.na(cur_hr) | cur_hr > desired_dim | is.na(cur_lr) | cur_lr < desired_dim) {
    message("Iteration: ", i, "\n", "Low Rank: ", lr[i], "\n High Rank: ",
            hr[i], "\n Target Rank: ",
            desired_dim)
    if(parallel == T) {
      if(is.null(high_gamma)) {
        l <- furrr::future_map(c(hk[i], lk[i]),
                               function(k)
                               {
                                 calc_gamma(Y, T_vec, J_vec, N, J, Time, Gamma = NULL, k = k,
                                            verbose = F)
                               })
      } else {
        l <- furrr::future_map2(list(hk[i], lk[i]),
                                list(high_gamma, low_gamma),
                                function(k, Gamma)
                                {
                                  calc_gamma(Y, T_vec, J_vec, N, J, Time, Gamma = Gamma, k = k,
                                             verbose = F)
                                })
      }

      high_gamma <- l[[1]]
      low_gamma <- l[[2]]

    } else {
      high_gamma <- calc_gamma(Y, T_vec, J_vec, N, J, Time, Gamma = high_gamma, k = hk[i],
                               verbose = verbose)
      low_gamma <- calc_gamma(Y, T_vec, J_vec, N, J, Time, Gamma = low_gamma, k = lk[i],
                              verbose = verbose)
    }


    hr[i] <- Matrix::rankMatrix(high_gamma)
    lr[i] <- Matrix::rankMatrix(low_gamma)

    cur_hr <- hr[i]
    cur_lr <- lr[i]

    # if rank under high k is too high (not <= to desired)
    # raise it because the desired rank is outside the bracket
    if(hr[i] > desired_dim) {
      hk[i + 1] <- hk[i] * 4
    } else {
      hk[i + 1] <- hk[i]
    }

    # if rank under low k is too low (not >= to desired)
    # lower it because the desired rank is outside the lower bracket
    if (lr[i] < desired_dim) {
      lk[i + 1] <- lk[i] / 4
    } else {
      lk[i + 1] <- lk[i]
    }
    i = i + 1
  }

  middle_gamma <- low_gamma

  while (stop == 0 & i <= maxiter) {

    if(verbose) {
      message("K Iteration: ", i, "\n",
              "Lower Penalization Parameter: ", signif(lk[i], 3), "\n",
              "Upper Penalization Parameter: ", signif(hk[i], 3), "\n",
              "Lower Penalization Rank: ", signif(lr[i-1], 3), "\n",
              "Upper Penalization Rank: ", signif(hr[i-1], 3), "\n")
    }
    # some version of the bisection method
    mk <- (hk[i] + lk[i])/2

    middle_gamma <- calc_gamma(Y, T_vec, J_vec, N, J, Time,
                               Gamma = middle_gamma, k = mk,
                               verbose = verbose)

    mr <- Matrix::rankMatrix(middle_gamma)

    rank_diff <- mr - desired_dim
    if(rank_diff > 0) {
      lk[i + 1] <- mk
      hk[i + 1] <- hk[i]
      low_gamma <- middle_gamma
      lr[i] <- Matrix::rankMatrix(low_gamma)
      hr[i] <- hr[i - 1]
    } else {
      hk[i + 1] <- mk
      lk[i + 1] <- lk[i]
      high_gamma <- middle_gamma
      hr[i] <- Matrix::rankMatrix(high_gamma)
      lr[i] <- lr[i - 1]
    }


    if(lr[i] == desired_dim) {
      return(low_gamma)
    }

    i = i + 1
  }

  Gamma
}

#' Generate predictions from OLS by group quickly
#' @param vec vector of estimated Gamma values, length of Y variable
#' @param groups vector of group variables to use in the independent regressions
#' @param reg_data matrix or data.frame of principal components, with one row for each gamma, one column for each desired latent factor
#' @param min_N minimum N per regression
#' @importFrom data.table data.table
#' @importFrom data.table `:=`
#' @import data.table
#' @export
ols_by_group <- function(vec, reg_data, groups, min_N){

  full_data <- data.table::data.table(vec = vec, Groups = groups, reg_data)
  full_data[,Order := .I]
  form <- paste0("Gamma ~ 0 + ", paste0("get(",colnames(reg_data),")", collapse = " + "))

  pred_ols <- function(l, y, W) {

    X = Reduce(cbind, l)
    if(length(y) <= ncol(X)){
      return(NA)
    } else
    as.numeric(X %*% solve(crossprod(X, X)) %*% crossprod(X, y))
  }

  return_data <- full_data[
        ,N := .N - sum(vec == 0), by = Groups
      ][N >= min_N, reconstructed_vec := pred_ols(l = .SD, y = vec), by = Groups,
        .SDcols = colnames(reg_data)
        ][N < min_N, reconstructed_vec := 0]

  data.table::setorder(return_data, Order)
  return_data[['reconstructed_vec']]

}

#' Iterate ols given Principal components regression data
#' @param e errors, ignored
#' @param est previous estimate, ignored
#' @param w weight vector
#' @param reg_data regression data from principal components
#' @param J_vec group vector to split ols
#' @param min_N minimum threshold for observations
gamma_ols_iter <- function(e, est, w, reg_data, J_vec, min_N, ...) {
  structure(ols_by_group(e * sqrt(w), reg_data,
                         J_vec, min_N), class = "group_ols")
}

#' "Predict" from group ols
#' @param original estimate
#' @param DT data.table, unused here
#' @param ... other parameters, ignored
predict.group_ols <- function(est, DT, ...) {
  est
}

#' Recover gamma through latent time and group factor structures
#' @param Y vector of outcomes
#' @param T_vec vector of time assignment
#' @param J_vec vector of group assignments
#' @param N number of assets
#' @param Time integer representing number of periods
#' @param J integer representing number of groups
#' @param tau quantile to calculate
#' @param rho stochastic relaxation parameter
#' @param err_u_thresh error threshold for u
#' @param min_N minimum N for group-wise OLS
#' @export
fe_calc_admm <- function(Y,
                         T_vec,
                         J_vec,
                         N,
                         J,
                         Time,
                         Gamma = NULL,
                         tau = 0.5,
                         rho = 0.1,
                         init_k = 1,
                         desired_dim = NULL,
                         err_u_thresh = 0.01,
                         err_g_thresh = 0.01,
                         verbose = F,
                         maxiter = 100,
                         min_N = 30) {



  T_vec <- factor(T_vec)
  T_levels <- levels(T_vec)

  if(is.null(Gamma)) {
    Gamma <- calc_reduce_rank_gamma(Y, T_vec, J_vec,
                                    N = N,
                                    Time = Time,
                                    J = J,
                                    verbose = verbose,
                                    tau = tau,
                                    init_k = init_k,
                                    desired_dim = desired_dim)
    Gamma <- as.matrix(Gamma)
  }


  pc_data <- data.table::as.data.table(stats::prcomp(t(Gamma))$x[,1:desired_dim])
  pc_data[, Year := factor(.I, labels = T_levels)]

  w <- data.table::data.table(Y, Group = J_vec)[,.(w = sqrt(.N)), by = Group]
  pc_cols <- setdiff(colnames(pc_data), "Year")
  reg_data <- data.table::data.table(Year = T_vec,
                                     Group = J_vec)[pc_data, on = "Year"][
                                       w, on = "Group"
                                     ]

  reg_data <- reg_data[,(pc_cols) := lapply(.SD, function(x) x * w),
                       .SDcols = pc_cols]

  admm(Y = Y, X = reg_data[,.(PC1, PC2, PC3)], est = rep(0, length(Y)),
       vec = rep(0, length(Y)),
       rho = rho, tau = tau,
       iterfun = gamma_ols_iter,
       k = k, err_u_thresh = err_u_thresh, verbose = verbose,
       reg_data = reg_data[,.(PC1, PC2, PC3)],
       J_vec = J_vec,
       min_N, w = reg_data$w)

}
