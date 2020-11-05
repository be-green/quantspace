#' Get coefficient Deltas to calculate Delta Method standard errors
#' @param spacing_coef coefficients from spacing regression
#' @param dep_col Column of response variable.
#' @param data Regression specification matrix.
#' @param alphas Quantiles to be estimated.
#' @param jstar First quantile to be estimated (usually the center one)
#' @param weights vector of observations weights. Defaults to 1/N
#'  for each observation if unspecified
#' @param desired_coefs which coefficients to calculate standard errors.
#'   If want standard errors for first and second parameters for each quantile,
#'   should desired_coeffs = c(1,2)
#' @param ... unused for now
#' @return List w/ p elements, with jth element a
#'        Nxk matrix corresponding to the delta_ij values
#'        for the jth quantile estimated
getDeltas <- function(spacing_coef, dep_col,data,alphas,jstar,weights=NULL,desired_coeffs=NULL, ...) {

  x <- data #Nxk
  y <- dep_col #Nx1
  quantiles <- spacingsToQuantiles(spacingCoef = spacing_coef,
                                   x, jstar) #Nxp

  N <- length(y) # Num observations
  k <- dim(x)[2] # Num parameters per quantile
  p <- length(alphas) # Num quantiles

  if (is.null(weights)){
    obswgt <- repMat(1/N,N,1)
  }
  if (!is.null(weights)){
    obswgt <- weights/sum(weights)
  }

  for (j in 1:p) {
    if (j == 1) {
      fyx <- eval_PDF(quantiles[,j],quantiles,alphas)
    }
    if (j > 1 ) {
      fyx <- cbind(fyx,eval_PDF(quantiles[,j],quantiles,alphas))
    }
  }

  #Initialize list of Delta matrices to an empty list of length p
  Delta <- rep(list(NaN),p)

  #For central quantile j*
  qjstar <- as.matrix(quantiles[,jstar]) #Nx1
  alphajstar <- alphas[jstar] #1x1
  mjstar <- x*(repMat((y < qjstar),1,k)-repMat(alphajstar,N,k)) #Nxk
  # mujstar <- repMat(t(sparse_Sums(as(repMat(obswgt,1,k)*mjstar,"dgCMatrix"))),N,1) #Nxk
  mujstar <- do.call(rbind,rep(list(spSums(repMat(obswgt,1,k)*mjstar)), N)) #Nxk

  Vjstar <- t(x)%*% spDiag(obswgt*fyx[,jstar]) %*% x #kxk

  Deltajstar <-  t(-solve(Vjstar)%*%t(mjstar - mujstar)) #Nxk

  #Add central quantile to Delta list in j* slot
  if (is.null(desired_coeffs)){
    Delta[[jstar]] <- Deltajstar
  }
  else {
    Delta[[jstar]] <- Deltajstar[,desired_coeffs]
  }

  #Values to start recursion
  Deltajold <- Deltajstar
  qjold <- qjstar
  Deltabar <- Deltajold
  epsold <- as.matrix(y - qjstar) #Nx1
  #Remove variables that are not needed and take up memory
  rm(qjstar,mjstar,mujstar,Vjstar)
  #Loop for upper quantiles j=j* + 1 to p
  if (jstar < p) {
    for (j in (jstar+1):p) {
      thetaj <-spacing_coef[,j] #kx1
      tauj <- (alphas[j] - alphas[j-1] )/(1-alphas[j-1]) #1x1
      mj <- x*(repMat((epsold < exp(x %*% thetaj)),1,k) - repMat(tauj,N,k))*repMat((epsold > 0),1,k) #Nxk
      muj <- do.call(rbind,rep(list(spSums(repMat(obswgt,1,k)*mj)), N))
      #muj <- repMat(t(spSums(repMat(obswgt,1,k)*mj)),N,1) #Nxk
      F2j <- fyx[,j] + fyx[,j-1]*repMat(tauj-1,N,1) #Nx1

      #Forming big Dj matrix
      for (jtilde in jstar:(j-1) ){
        thetajtilde <-spacing_coef[,jtilde] #kx1
        if ( jtilde == jstar){
          Dj <- t(x) %*% spDiag(obswgt*F2j) %*% x #kxk
        }
        if (jtilde > jstar) {
          D <- t(x) %*% spDiag(obswgt*F2j*exp(x %*% thetajtilde)) %*% x   #kxk
          Dj <- rbind(Dj,D)
        }
      }
      Dj <- t(Dj) #kx(j-jstar)*k

      fyx_expxtheta_eps <- fyx[,j]*exp(x %*% thetaj) #Nx1
      Vj <- t(x) %*% spDiag(obswgt*fyx_expxtheta_eps) %*% x #kxk
      Deltaj <-  t(-solve(Vj)%*%(t(mj - muj) + Dj%*%t(Deltabar))) #Nxk

      #Values to use in next loop of recursion
      Deltabar <- cbind(Deltabar,Deltaj) #(j+1-jstar)*kxN
      qjold <- quantiles[,j]
      epsold <- as.matrix(y - qjold)
      Deltajold <- Deltaj

      #Add Deltaj to Delta list in jth slot
      #Add Deltaj to Delta list in jth slot
      if (is.null(desired_coeffs)){
        Delta[[j]] <- Deltaj
      }
      else {
        Delta[[j]] <- Deltaj[,desired_coeffs]
      }
      #Remove variables that are not needed and take up memory
      # rm(thetaj,mj,muj,F2j,thetajtilde,fyx_expxtheta_eps,Vj,Deltaj)
    }
  }

  #Loop for lower quantiles j*-1 to 1
  #Values to start recursion
  Deltajold <- Deltajstar
  qjold <- quantiles[[jstar]]
  Deltabar <- Deltajold
  epsold <- as.matrix(y - quantiles[[jstar]]) #Nx1
  rm(Deltajstar)
  if (jstar > 1) {
    for (j in (jstar-1):1){
      thetaj <-spacing_coef[,j] #kx1
      tauj <- (alphas[j+1] - alphas[j] )/(alphas[j+1]) #1x1
      mj <- x*(repMat((epsold > -exp(x %*% thetaj)),1,k) -
                 repMat(tauj,N,k))*repMat((epsold < 0),1,k) #Nxk
      muj <- do.call(rbind,rep(list(spSums(repMat(obswgt,1,k)*mj)), N))
      F2j <- fyx[,j] + fyx[,j+1]*repMat(tauj-1,N,1) #Nx1

      #Forming big Dj matrix
      for (jtilde in jstar:(j+1) ){
        thetajtilde <-spacing_coef[,jtilde] #kx1
        if ( jtilde == jstar){
          Dj <- -t(x) %*% spDiag(obswgt*F2j) %*% x #kxk
        }
        if (jtilde < jstar) {
          D <- t(x) %*% spDiag(obswgt*F2j*exp(x %*% thetajtilde)) %*% x   #kxk
          Dj <- rbind(Dj,D)
        }
      }
      Dj <- t(Dj) #kx(j-jstar)*k

      fyx_expxtheta_eps <- fyx[,j]*exp(x %*% thetaj) #Nx1
      Vj <- t(x) %*% spDiag(obswgt*fyx_expxtheta_eps) %*% x #kxk
      Deltaj <-  t(-solve(Vj)%*%(t(mj - muj) + Dj%*%t(Deltabar))) #Nxk

      #Values to use in next loop of recursion
      Deltabar <- cbind(Deltabar,Deltaj) #(jstar-(j+1))*kxN
      qjold <- quantiles[,j]
      epsold <- as.matrix(y - qjold)
      Deltajold <- Deltaj

      #Add Deltaj to Delta list in jth slot
      if (is.null(desired_coeffs)){
        Delta[[j]] <- Deltaj
      }
      else {
        Delta[[j]] <- Deltaj[,desired_coeffs]
      }
      #Remove variables that are not needed and take up memory
      # rm(thetaj,mj,muj,F2j,thetajtilde,fyx_expxtheta_eps,Vj,Deltaj)
    }
  }

  Delta
}

#' Computes variance-covariance matrix and standard errors of parameter estimates
#' @param delta list of length p, where each element is
#' an Nxk matrix of delta_ij values with jth element corresponding to jth quantile
#' estimated (would be obtained by running getSEDeltas function)
#' @param cluster_indices vector of length N corresponding to labels
#' of which cluster given observation belongs to (defaults to NULL)
#' weights--observation weights. Defaults to NULL, which corresponds
#' to using 1/N as the observation weights
#' @param ... unused for now
#' @return A list with two elements: varcov--kp X kp variance covariance matrix of parameters,
#' se--kp x 1 standard errors of parameters
getDeltaVarCov <- function(delta,cluster_indices=NULL,weights=NULL, ...){

  #Combine list of deltas into one large matrix
  for (i in 1:length(delta)){
    if (i == 1) {
      deltaCombined <- delta[[i]]
    }
    if (i > 1){
      deltaCombined <- cbind(deltaCombined,delta[[i]])
    }
  }
  N = dim(deltaCombined)[1] #Number of observations
  #Case where weights are not provided
  if (is.null(weights)){
    obswgt <- repMat(1/N,N,1)
  }
  #Case where weights are provided
  if (!is.null(weights)){
    obswgt <- as.matrix(weights/sum(weights))
  }

  #Case where no cluster indices are specified
  if (is.null(cluster_indices)){
    var_est = t(deltaCombined) %*% spDiag(obswgt) %*% deltaCombined
  }

  #Case with clustering
  if (!is.null(cluster_indices)) {

    nparams <- dim(deltaCombined)[2]
    nclust <- dim(cluster_indices)[1]

    #This is much quicker and better
    for (i in 1:nclust){
      if (i==1){
        clustindex <- i*(matrix(1,1+cluster_indices[i,2]-cluster_indices[i,1],1))
      }
      else{
        clustindex <- rbind(clustindex,
                            i*(matrix(1,1+cluster_indices[i,2]-cluster_indices[i,1],1)) )
      }
    }

    DATA <- as(repMat(sqrt(obswgt),1,nparams)*deltaCombined,"dgCMatrix")
    D <- aggregate.Matrix(DATA,groupings = as.factor(clustindex),FUN=spSums)
    var_est <- t(D)%*%D
  }

  #Return estimated variance-covariance matrix of parameters
  #and standard error estimates
  return(list('varcov'=var_est,'se'=sqrt(diag(var_est)/N)))
}

