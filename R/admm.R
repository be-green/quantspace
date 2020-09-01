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

#' Calculate Gamma
#' @param DT data.table of values to be used in optimization
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

  if(is.null(Gamma)) {
    Gamma <- data.table::as.data.table(matrix(0, J, Time))
  }

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

  # Big matrix of ones if needed
  One <- matrix(1, NT, 1)

  err_u <- 1000
  iter <- 1
  DT[,Gamma := 0]
  nu <- k*max(sqrt(N), sqrt(Time))


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

  while(err_u > err_u_thresh) {

    if (iter == 1) {
      # Initial u
      DT[,u := rep(0, NT)]
    } else {

      err_u <- norm(as.matrix(rho * (DT$Y - DT$r - DT$Gamma)))/(sum(abs(DT$u)) + 1)
      DT[,u := u + rho * (Y - r - Gamma)]

      if(verbose) {
        message(paste0("Iteration: ", iter, "\n",
                       "Error: ", signif(err_u, 4)))
      }

    }

    # Append to data frame
    DT[,r := calc_r(rho, u, Y, Gamma, tau, One)]
    DT[,e := Y - r + u / rho]

    # Ugly but fast
    E <- dcast(
      DT[,.(e1 = sum(e)/length(e)),
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

    merge_gamma <- data.table::melt(data.table::copy(Gamma)[,Group := factor(1:.N, labels = levels(J_levels))],
                                    id.vars = "Group", variable.factor = F, variable.name = "Year",
                                    value.name = "Gamma")[,
                                                          Year := factor(as.integer(Year),
                                                                         labels = T_levels)]

    DT <- merge_gamma[DT, on = c("Group", "Year")]
    iter = iter + 1
  }
  return(t(Gamma))
}

