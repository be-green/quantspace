quad_form <- function(a, b, c, discriminant = NULL) {
  if (!is.null(discriminant)) {
    root1 <- (-b + sqrt(discriminant)) / (2*a)
    root2 <- (-b - sqrt(discriminant)) / (2*a)
  } else {
    root1 <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    root2 <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
  }
  return(cbind(root1, root2))
}

cub_root <- function(a, b, c, d) {
  roots <- apply(cbind(d, c, b, a), 1, polyroot)
  return(Re(roots)[abs(Im(roots)) < 1e-6])
}

cub_root_select <- function(a, b, c, d, q_1, q_2) {
  roots <- apply(cbind(d, c, b, a), 1, polyroot)
  roots_clean <- matrix(Re(roots)[abs(Im(roots)) < 1e-6], ncol = 3, byrow = TRUE)
  select <- mapply(function(x, q_1, q_2) q_1 <= x && x <= q_2, roots_clean, q_1, q_2)
  return(roots_clean[select])
}

cub_root_select_rconics <- function(a, b, c, d, q_1, q_2) {
  #roots_clean <- matrix(Re(roots)[abs(Im(roots)) < 1e-6], ncol = 3, byrow = TRUE)
  coeff <- cbind(a,b,c,d)
  roots <- apply(cbind(a,b,c,d), 1, cubic)
  roots_clean <- matrix(roots, ncol = 3, byrow = TRUE)
  roots_clean[abs(Im(roots_clean)) > 1e-6] <- NA
  roots_clean <- Re(roots_clean)
  #select <- mapply(function(x, q_1, q_2) q_1-10e-6 <= x && x <= q_2+10e-6, roots_clean, q_1, q_2)
  select <- mapply(function(x, q_1, q_2) q_1-1e-6 <= x && x <= q_2+1e-6, roots_clean, q_1, q_2)
  roots_clean[!select] <- NA
  roots_clean <- apply(roots_clean, 1, min, na.rm = TRUE)
  return(roots_clean)
}

cub_root_deriv <- function(a, b, c, d) {
  roots <- apply(cbind(d, c, b, a), 1, polyroot)
  root <- Re(roots)[abs(Im(roots)) < 1e-6]
  return(1 / (c + 2 * root * b + 3 * root^2 * a))
}

q_spline_R <- function(quantiles, alphas, tails = "gaussian") {
  ###########################################################
  #A function to compute the tail parameters and the second derivatives
  #Input:
  #quantiles: a matrix of size N b p;  the fitted quantiles computed
  #           based alphas
  #tails: based on what distribution to compute the tails
  #Retrun:
  #       A list containing the tail parameters and the second
  #       derivatives
  ###########################################################

  p <- dim(quantiles)[2]
  N <- dim(quantiles)[1]

  q_shift <- quantiles[, 1]
  q_stretch <- quantiles[, p] - quantiles[,1]
  #q_shift <- rep(0,1000)
  #q_stretch <- rep(1,1000)
  quantiles <- (quantiles - q_shift) / q_stretch

  # solving for location and scale parameters for a distribtuion for the left and right tails
  if (tails == "gaussian") {
    tail_param_l <- quantiles[, 1:2    ] %*% t(inv(matrix(c(1, 1, qnorm(alphas[1  ]), qnorm(alphas[2])), nrow = 2, ncol = 2)))
    tail_param_u <- quantiles[, (p-1):p] %*% t(inv(matrix(c(1, 1, qnorm(alphas[p-1]), qnorm(alphas[p])), nrow = 2, ncol = 2)))

    # calculating derivatives of the pdf -- boundary conditions for the spline
    yp1 <- dnorm(quantiles[, 2  ], mean = tail_param_l[,1], sd = tail_param_l[,2])
    ypp <- dnorm(quantiles[, p-1], mean = tail_param_u[,1], sd = tail_param_u[,2])
  } else if (tails == "exponential") {
    tail_param_l <- quantiles[, 1:2    ] %*% t(inv(matrix(c(1, 1, log(alphas[1]), log(alphas[2])), nrow = 2, ncol = 2)))
    tail_param_u <- quantiles[, (p-1):p] %*% t(inv(matrix(c(1, 1, -log(1 - alphas[p-1]), -log(1 - alphas[p])), nrow = 2, ncol = 2)))

    # calculating derivatives of the pdf -- boundary conditions for the spline
    yp1 <- exp(1)^( (quantiles[, 2]   - tail_param_l[, 1])/tail_param_l[, 2]) / tail_param_l[, 2]
    ypp <- exp(1)^(-(quantiles[, p-1] - tail_param_u[, 1])/tail_param_u[, 2]) / tail_param_u[, 2]
  } else {
    stop("Passed value to tails argument must be either 'gaussian' or 'exponential'.")
  }
  # print(matrix(c(1, 1, qnorm(alphas[1  ]), qnorm(alphas[2])), nrow = 2, ncol = 2))

  # next, we solve for the parameters of the cubic spline characterizing the CDF. Per the numerical recipes text,
  # this can be characterized in terms of the second derivative at each knot

  # from this point on, the interpolation script doesn't need the lowest and highest quantiles. Dropping for
  # more intuitive indexing
  quantiles <- matrix(quantiles[, -c(1,p)], nrow = N)
  alphas <- alphas[-c(1,p)]
  p <- p - 2

  # this is the tridiagonal algorithm from the numerical recipes book
  if(p>1){
    u <- matrix(rep(1, N*(p-1)), nrow = N, ncol = p-1)
    y2 <- matrix(nrow = N, ncol = p)

    y2[,1] <- -0.5
    u[,1] <- (3 / (quantiles[,2] - quantiles[,1])) * ((alphas[2]-alphas[1])/(quantiles[,2]-quantiles[,1]) - yp1)
    for(i in 2:(p-1)) {
      sig <- (quantiles[,i]-quantiles[,i-1]) / (quantiles[,i+1]-quantiles[,i-1])
      d <- sig*y2[,i-1] + 2
      y2[,i] <- (sig-1) / d
      u[,i] <- (alphas[i+1]-alphas[i]) / (quantiles[,i+1]-quantiles[,i]) - (alphas[i]-alphas[i-1]) / (quantiles[,i]-quantiles[,i-1])
      u[,i] <- (6*u[,i]/(quantiles[,i+1]-quantiles[,i-1]) - sig*u[,i-1]) / d
    }
    qp <- 0.5
    up <- (3 / (quantiles[,p]-quantiles[,p-1])) * (ypp - ((alphas[p]-alphas[p-1])/(quantiles[,p]-quantiles[,p-1])))

    y2[,p] <- (up - (qp*u[,p-1])) / ((qp*y2[,p-1]) + 1)
    for(i in (p-1):1) {
      y2[,i] <- y2[,i]*y2[,i+1]+u[,i]
    }
  }
  else{
    y2<-NULL}

  return(list(y2 = y2, tail_param_u = tail_param_u, tail_param_l = tail_param_l, q_shift = q_shift, q_stretch = q_stretch))

}

splint_R <- function(y, quantiles, alphas, y2, tail_param_u, tail_param_l,
                     q_shift, q_stretch, y_rows = NULL,
                     tails = "gaussian", distn = "c") {
  ###########################################################
  #A function to condut the interpolation given data and fitted quantiles
  #y: a vector of quantile or pdf
  #quantiles: a matrix of size N b p;  the fitted quantiles computed
  #           based alphas
  #y2: second derivatives
  #tail_param_u: the tail parameters for upper quantiles
  #tail_param_l: the tail parameters for lower quantiles
  #q_shift: normalization parameters
  #q_stretch: normalization parameters
  #tails: based on what distribution to compute the tails
  #distn: a string: takes on value:
  #       p: pdf
  #       c: cdf
  #       q: quantile
  #Retrun:
  #       y_hat: A vector of quantiles or density corresponding to y
  ###########################################################
  N <- length(y)
  p <- length(alphas)

  # applying location transformations defined earlier
  quantiles <- (quantiles - q_shift) / q_stretch

  # as before, we don't need the extreme quantiles for anything. Dropping for easier indexing
  quantiles <- matrix(quantiles[, -c(1,p)], nrow = N)
  alphas <- alphas[-c(1,p)]
  p <- p - 2


  # Make the distribution argument into a vector
  if (length(distn) == 1) {
    distn <- rep(distn, N)
  }

  # Apply location and scale normalizations
  y <- ifelse(distn != "q",
              (y - q_shift) / q_stretch,
              y)

  ####### These variables can probably be cleaned up a bit ########
  y_hat_prime <- vector(mode = "numeric", length = N)

  # preallocating vectors here
  y_hat <- numeric(N)
  y_hat_low <- rep(NA, N)
  y_hat_hi <- rep(NA, N)
  y_hat_mid <- rep(NA, N)

  #setting tails to missing where
  y_low <- y
  y_hi <- y
  y_low[y > alphas[1] & distn == "q"] <- NA
  y_hi[y < alphas[p] & distn == "q"] <- NA
  y_low[y > quantiles[, 1] & distn != "q"] <- NA
  y_hi[y < quantiles[, p] & distn != "q"] <- NA

  # If the observation is in the tails, we just need to evaluate the closed form CDF/PDF/quantile function
  if (tails == "gaussian") {
    if (length(distn[!is.na(y_low)]) == 0) {
      y_hat_low <- rep(NA, N)
    } else {
      y_hat_low <- ifelse(distn == "c",
                          pnorm(y_low, mean = tail_param_l[, 1], sd = tail_param_l[, 2]),
                          ifelse(distn == "p",
                                 dnorm(y_low, mean = tail_param_l[, 1], sd = tail_param_l[, 2]),
                                 qnorm(y_low, mean = tail_param_l[, 1], sd = tail_param_l[, 2])
                          )
      )
    }
    if (length(distn[!is.na(y_hi)]) == 0) {
      y_hat_hi <- rep(NA, N)
    } else {
      y_hat_hi <- ifelse(distn == "c",
                         pnorm(y_hi, mean = tail_param_u[, 1], sd = tail_param_u[, 2]),
                         ifelse(distn == "p",
                                dnorm(y_hi, mean = tail_param_u[, 1], sd = tail_param_u[, 2]),
                                qnorm(y_hi, mean = tail_param_u[, 1], sd = tail_param_u[, 2])
                         )
      )
    }

  } else if (tails == "exponential") {
    if (length(distn[!is.na(y_low)]) == 0) {
      y_hat_low <- rep(NA, N)
    } else {
      y_hat_low <- ifelse(distn == "c",
                          exp(1)^( (y_low - tail_param_l[, 1])/tail_param_l[, 2]),
                          exp(1)^( (y_low - tail_param_l[, 1])/tail_param_l[, 2]) / tail_param_l[, 2]
      )
    }
    if (length(distn[!is.na(y_hi)]) == 0) {
      y_hat_hi <- rep(NA, N)
    } else {
      y_hat_hi <- ifelse(distn == "c",
                         1-exp(1)^(-(y_hi  - tail_param_u[, 1])/tail_param_u[, 2]),
                         exp(1)^(-(y_hi  - tail_param_u[, 1])/tail_param_u[, 2]) / tail_param_u[, 2]
      )
    }
  } else {
    stop("Passed value to tails argument must be either 'gaussian' or 'exponential'.")
  }

  # next, we execute a block to do the interior of the spline
  if(p>1){
    y_mid <- y
    y_mid[(y <= quantiles[, 1] | y >= quantiles[, p]) & distn != "q"] <- NA
    y_mid[(y <= alphas[1] | y >= alphas[p]) & distn == "q"] <- NA

    for (i in 1:(p-1)) {
      y_int <- y_mid >= quantiles[, i] & y_mid <= quantiles[, i+1] & distn != "q" | y_mid >= alphas[i] & y_mid <= alphas[i+1] & distn == "q"
      y_int[is.na(y_int)] <- FALSE
      if (!any(y_int)) {
        next
      }
      # print(i)
      h <- quantiles[y_int, i+1] - quantiles[y_int, i]
      a <- (quantiles[y_int, i+1] - y_mid[y_int])/h
      b <- (y_mid[y_int] - quantiles[y_int, i])/h
      b_prime <- 1/h
      a_prime <- -b_prime
      y_hat_mid_temp <- vector(mode = "numeric", length = length(h))
      discriminant <- (1/(3*h^2))*(-6*alphas[i]*y2[y_int, i] + 6*alphas[i+1]*y2[y_int, i] + h^2*y2[y_int, i]^2 + 6*alphas[i]*y2[y_int, i+1] - 6*alphas[i+1]*y2[y_int, i+1] - 2*h^2*y2[y_int, i]*y2[y_int, i+1] + 3*quantiles[y_int, i]^2*y2[y_int, i]*y2[y_int, i+1] - 6*quantiles[y_int, i]*quantiles[y_int, i+1]*y2[y_int, i]*y2[y_int, i+1] + 3*quantiles[y_int, i+1]^2*y2[y_int, i]*y2[y_int, i+1] + h^2*y2[y_int, i+1]^2)
      zeros2 <- ifelse(cbind(discriminant,discriminant) >= 0,
                       quad_form(-(y2[y_int, i]/(2*h)) + y2[y_int, i+1]/(2*h), (quantiles[y_int, i+1]*y2[y_int, i])/h - (quantiles[y_int, i]*y2[y_int, i+1])/h, -(alphas[i]/h) + alphas[i+1]/h + (h*y2[y_int, i])/6 - (quantiles[y_int, i+1]^2*y2[y_int, i])/(2*h) - (h*y2[y_int, i+1])/6 + (quantiles[y_int, i]^2*y2[y_int, i+1])/(2*h), discriminant = discriminant),
                       NA
      )
      z_in_int <- zeros2 >= quantiles[y_int, i] & zeros2 <= quantiles[y_int, i+1]

      #note that z_in_int_comb is TRUE if there is a problem with interpolation
      z_in_int_comb <- ifelse(is.na(zeros2[,1]), FALSE, z_in_int[,1] | z_in_int[,2])


      a1 <- (y2[y_int, i+1] - y2[y_int, i]) / (6 * h);
      b1 <- (quantiles[y_int, i+1] * y2[y_int, i] - quantiles[y_int, i] * y2[y_int, i+1]) / (2 * h);
      c1 <- (alphas[i+1] - alphas[i]) / h + (h * (y2[y_int, i] - y2[y_int, i+1])) / 6 + (quantiles[y_int, i] * quantiles[y_int, i] * y2[y_int, i+1] - quantiles[y_int, i+1] * quantiles[y_int, i+1] * y2[y_int, i]) / (2 * h);
      d1 <- (alphas[i] * quantiles[y_int, i+1] - alphas[i+1] * quantiles[y_int, i]) / h + (quantiles[y_int, i+1] * quantiles[y_int, i+1] * quantiles[y_int, i+1] * y2[y_int, i] - quantiles[y_int, i] * quantiles[y_int, i] * quantiles[y_int, i] * y2[y_int, i+1]) / (6 * h) + (h * (quantiles[y_int, i] * y2[y_int, i+1] - quantiles[y_int, i+1] * y2[y_int, i])) / 6;

      # evaluate expressions here. If the spline is not monotonic, z_in_int_comb is true. Then, we use a linear interpolation instead
      # otherwise, we use the spline

      y_hat_mid_temp <- ifelse(z_in_int_comb, {
        ifelse(distn[y_int] == "c",
               (y_mid[y_int] - quantiles[y_int, i]) * (alphas[i+1] - alphas[i]) / h + alphas[i],
               ifelse(distn[y_int] == "p",
                      (alphas[i+1] - alphas[i]) / h,
                      (y_mid[y_int] - alphas[i]) * h / (alphas[i+1] - alphas[i]) + quantiles[y_int, i]
               )
        )
      }, {
        ifelse(distn[y_int] == "c",
               a*alphas[i] + b*alphas[i+1] + ((a^3 - a)*y2[y_int, i] + (b^3 - b)*y2[y_int, i+1]) * (h^2) / 6,
               ifelse(distn[y_int] == "p",
                      a_prime*alphas[i] + b_prime*alphas[i+1] + ((3*(a^2)*a_prime - a_prime)*y2[y_int, i] + (3*(b^2)*b_prime - b_prime)*y2[y_int, i+1]) * (h^2) / 6,
                      cub_root_select_rconics(a1, b1, c1, d1 - y_mid[y_int], quantiles[y_int, i], quantiles[y_int, i+1])
               )
        )
      })
      y_hat_mid[y_int] <- y_hat_mid_temp
    }
  }
  # combining the three different calculations here
  y_hat <- pmin(y_hat_mid, y_hat_low, na.rm = TRUE)
  y_hat <- pmin(y_hat, y_hat_hi, na.rm = TRUE)

  # Undo location scal normalizations from earlier
  y_hat <- ifelse(distn == "q",
                  y_hat * q_stretch + q_shift,
                  ifelse(distn == "p",
                         y_hat / q_stretch,
                         y_hat)
  )

  return(y_hat)
}


splint <- function(y, quantiles, alphas, y2, tail_param_u, tail_param_l, y_rows = NULL, tails = "gaussian", distn = "c") {
  if (distn == "c") {
    splint_R(y, quantiles, alphas, y2, tail_param_u, tail_param_l, y_rows, tails, TRUE)
  } else {
    splint_R(y, quantiles, alphas, y2, tail_param_u, tail_param_l, y_rows, tails, FALSE)

  }
}


eval_density_R<-function(y, alphas, m = 1, s = 1, distn = "p"){
  ###########################################################
  #A function to evaluate the quantile or density of
  #given data based on normal distribution
  #y: a vector of quantile or pdf
  #alphas: a vector of specifed quantiles to interpolate
  #m: mean of the normal distribution
  #s: standard deviation of the normal distribution
  #distn: a string: takes on value:
  #       p: pdf
  #       c: cdf
  #       q: quantile
  #Return:
  #       A vector of quantiles or density corresponding to y
  ###########################################################
  p<-length(alphas)
  jstar<-round(p/2)
  N<-length(y)

  #####set up the quantile parameters
  eta_theta_temp <- log(qnorm(alphas[-1], mean = m, sd = s) - qnorm(alphas[-length(alphas)], mean = m, sd = s))
  eta_theta <- matrix(c(eta_theta_temp[1:(round(p/2)-1)], m, eta_theta_temp[round(p/2):(p-1)]), nrow = 1)

  quantiles<-spacingsToQuantiles(eta_theta, matrix(1, N, 1), jstar)
  spline_params <- q_spline_R(quantiles, alphas, tails = "gaussian")
  y_hat <- splint_R(y, quantiles, alphas, spline_params[["y2"]], spline_params[["tail_param_u"]], spline_params[["tail_param_l"]],  spline_params[["q_shift"]], spline_params[["q_stretch"]],tails = "gaussian", distn = distn)
  return(y_hat)
}

eval_PDF <- function(y,quantiles,alphas,distn = "p", tails = "gaussian") {
  spline_params <- q_spline_R(quantiles,alphas,tails)
  y_hat <- splint_R(y, quantiles, alphas, spline_params[["y2"]], spline_params[["tail_param_u"]], spline_params[["tail_param_l"]],  spline_params[["q_shift"]], spline_params[["q_stretch"]],tails = tails, distn = distn)
  return(y_hat)
}

eval_CDF <- function(y,quantiles,alphas,distn = "c", tails = "gaussian") {
  spline_params <- q_spline_R(quantiles,alphas,tails)
  y_hat <- splint_R(y, quantiles, alphas, spline_params[["y2"]], spline_params[["tail_param_u"]], spline_params[["tail_param_l"]],  spline_params[["q_shift"]], spline_params[["q_stretch"]],tails = tails, distn = distn)
  return(y_hat)
}

eval_Quantiles <- function(y,quantiles,alphas,distn = "q", tails = "gaussian") {
  spline_params <- q_spline_R(quantiles,alphas,tails)
  y_hat <- splint_R(y, quantiles, alphas, spline_params[["y2"]], spline_params[["tail_param_u"]], spline_params[["tail_param_l"]],  spline_params[["q_shift"]], spline_params[["q_stretch"]],tails = tails, distn = distn)
  return(y_hat)
}
