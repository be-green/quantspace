// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::plugins(cpp11)]]

double checkfun (double x, double tau) {
  return x * (tau - (x < 0));
}

double checkfun (arma::vec& res, double tau, int n) {
  double value = 0;
  for(int i = 0; i < n; i++) {
    value += checkfun(res(i), tau);
  }
  return value;
}

// draw random exponential weights for bootstrap
// with rate 1
// [[Rcpp::export]]
arma::vec fast_rexp(int n) {
  // quantile function for rate 1
  return -1 * arma::log(1 - arma::randu(n));
}

// [[Rcpp::export]]
void update_huber_grad(const arma::mat& X_t,
                       const arma::vec& res,
                       arma::vec& derivs,
                       arma::vec& grad,
                       double tau,
                       double mu,
                       int n,
                       double one_over_n) {
  for(int i = 0; i < n; i = i + 1) {
    if (res(i) > mu) {
      derivs(i) = tau;
    } else if(res(i) < -mu) {
      derivs(i) = tau - 1;
    } else if(res(i) >= 0) {
      derivs(i) = res(i) * tau / mu;
    } else {
      derivs(i) = res(i) * (tau - 1) / mu;
    }
  }
  grad = X_t * derivs * one_over_n;
}

// [[Rcpp::export]]
void z_score(arma::mat& X, const arma::rowvec& mx, const arma::vec& sx, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) / sx(i);
  }
}

void reorder_columns(arma::mat& X, int intercept) {
  if(intercept != 0) {
    intercept = intercept - 1;
    if(intercept > 0) {
      X.insert_cols(intercept, 1);
      X.shed_col(0);
      X.col(intercept) = arma::ones(X.n_rows);
    }
  }
}

// [[Rcpp::export]]
arma::vec huber_grad_descent(const arma::colvec& y, const arma::mat& X,
                             const arma::mat& X_t, arma::vec& beta,
                             double tau, double n, double one_over_n,
                             int p, int maxiter, double mu,
                             double beta_tol, double check_tol) {


  // gradient vector
  arma::vec grad(p);
  arma::vec last_beta = beta;
  arma::vec last_grad = grad;
  arma::vec derivs(n);

  // vector of residuals
  arma::vec resid = y - X * beta;

  // differences from previous gradients
  // and betas
  arma::vec beta_diff(p);
  arma::vec grad_diff(p);

  double checkfun_diff = checkfun(resid, tau, n);
  double last_checkfun = checkfun_diff;
  double this_checkfun;

  double i = 1;
  double cross = 0;
  double delta = std::min(1/tau, 1/(1 - tau));
  // inf is the "max" norm over the vector
  while((i < maxiter) && ((arma::norm(grad, "inf") > beta_tol) || (i == 1)) &&
        (checkfun_diff / delta > check_tol || delta < 0.1)) {


    // step size and original gradient descent code from conquer package
    // https://github.com/XiaoouPan/conquer
    // picking step size
    delta = 1;
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(grad_diff.t() * grad_diff);
      double a2 = arma::as_scalar(beta_diff.t() * beta_diff) / cross;

      // pick the smaller and don't go _really_ off the rails
      delta = std::min(std::min(a1, a2), 1.0);
    }

    if(std::fmod(i, 100) == 0) {
      Rcpp::checkUserInterrupt();
    }

    // store last gradient for step size selection
    last_grad = grad;
    update_huber_grad(X_t, resid, derivs, grad, tau, mu, n, one_over_n);

    beta_diff = beta - last_beta;
    grad_diff = grad - last_grad;

    beta += (i - 1)/(i + 2) * beta_diff;
    beta += delta * grad;

    // update the beta movement for step size later
    beta_diff += delta * grad;

    // update residual vector
    resid -= X * beta_diff;
    last_beta = beta;

    this_checkfun = checkfun(resid, tau, n);
    checkfun_diff = abs((last_checkfun - this_checkfun));
    last_checkfun = this_checkfun;
    i++;
    cross = arma::as_scalar(beta_diff.t() * grad_diff);
  }
  return beta;
}


//' Compute quantile regression via accelerated gradient descent using
//' Huber approximation, warm start based on data subset
//' @param y outcome vector
//' @param X design matrix
//' @param X_sub subset of X matrix to use for "warm start" regression
//' @param y_sub subset of y to use for "warm start" regression
//' @param tau target quantile
//' @param init_beta initial guess at beta
//' @param intercept location of the intercept column, using R's indexing
//' @param num_samples number of samples used for subset of matrix used for warm start
//' @param mu neighborhood over which to smooth
//' @param maxiter maximum number of iterations to run
//' @param check_tol loss function change tolerance for early stopping
//' @param beta_tol tolerance for largest element of gradient, used
//' for early stopping
//' @param warm_start integer indicating whether to "warm up" on a subsample
//' of the data
//' @param scale whether to scale x & y variables
//' @export
// [[Rcpp::export]]
arma::vec fit_approx_quantile_model(arma::mat& X,
                                    arma::vec& y,
                                    arma::mat& X_sub,
                                    arma::vec& y_sub,
                                    double tau,
                                    arma::colvec init_beta,
                                    double mu = 1e-15,
                                    int maxiter = 100000,
                                    double beta_tol = 1e-4,
                                    double check_tol = 1e-6,
                                    int intercept = 1,
                                    double num_samples = 1000,
                                    int warm_start = 1,
                                    int scale = 1) {

  // p is dim, n is obs, one_over_n to avoid repeated calcs
  int p = X.n_cols;
  double n = X.n_rows;

  arma::vec q = {tau};
  // calc'd here because we use it a bunch
  double one_over_n = 1/n;

  arma::rowvec mx;
  arma::vec sx;
  double my;
  double sy;

  if(scale == 1) {

    // standardizing starts
    if(intercept > 0) {
      // remove intercept, we will calculate after
      // the subtraction is just to match R's indexing
      X.shed_col(intercept - 1);
      if(warm_start == 1) {
        X_sub.shed_col(intercept - 1);
      }
    }

    // standardizing everything to work w/ z scores
    // we will transform betas back at the end
    mx = arma::mean(X, 0);
    sx = arma::stddev(X, 0, 0).t();

    // standardize X
    for (int i = 0; i < X.n_cols; i++) {
      X.col(i) = (X.col(i) - mx(i)) / sx(i);

      if(warm_start == 1) {
        X_sub.col(i) = (X_sub.col(i) - mx(i)) / sx(i);
      }
    }

    my = arma::mean(y);
    sy = arma::stddev(y);
    if(intercept > 0) {
      y -= my;
      if(warm_start == 1) {
        y_sub -= my;
      }
    }
  }

  // put intercept at beginning after z-scoring
  // will re-arrange to match design matrix
  // at the end of the function
  if(intercept > 0) {
    X = arma::join_rows(arma::ones(n), X);
    if(warm_start == 1) {
      X_sub = arma::join_rows(arma::ones(num_samples),X_sub);
    }
  }
  // pre-transposed X
  arma::mat X_t = arma::trans(X);


  if(warm_start == 1) {

    // warm start
    arma::mat X_t_sub = arma::trans(X_sub);

    double one_over_num_samples = 1/num_samples;

    init_beta = huber_grad_descent(y_sub,
                                   X_sub,
                                   X_t_sub,
                                   init_beta,
                                   tau, num_samples,
                                   one_over_num_samples, p,
                                   100, mu,
                                   beta_tol, check_tol);

  }

  init_beta(0) = arma::as_scalar(arma::quantile(y - X.cols(1, p - 1) * init_beta.rows(1, p - 1),
                                 q));

  // full data
  arma::vec beta = huber_grad_descent(y, X, X_t, init_beta,
                                      tau, n, one_over_n, p,
                                      maxiter, mu,
                                      beta_tol, check_tol);

  if(scale == 1) {

    // unstandardize
    y += my;
    double m;
    double s;

    // unstandardize coefficients
    if(intercept > 0) {
      for (int i = 1; i < X.n_cols; i++) {
        m = mx(i - 1);
        s = sx(i - 1);
        X.col(i) = (X.col(i) + m) * s;
      }
      beta.rows(1, p - 1) /= sx;
      beta(0) += my - arma::as_scalar(mx * beta.rows(1, p - 1));
    } else {

      for (int i = 0; i < p; i++) {
        m = mx(i);
        s = sx(i);
        X.col(i) = (X.col(i) + m) * s;
      }

      beta.rows(0, p - 1) /= sx;
    }
    // if the intercept column of X isn't the first column,
    // re-order the coefficients
    if(intercept > 1) {
      reorder_columns(X, intercept);
      beta.insert_rows(intercept - 1, beta(0));
      beta.shed_row(0);
    }
  }
  return(beta);
}

void update_huber_grad_sp(const arma::sp_mat& X_t,
                         const arma::vec& res,
                         arma::vec& derivs,
                         arma::vec& grad,
                         double tau,
                         double mu,
                         int n,
                         double one_over_n) {
  for(int i = 0; i < n; i = i + 1) {
    if (res(i) > mu) {
      derivs(i) = tau;
    } else if(res(i) < -mu) {
      derivs(i) = tau - 1;
    } else if(res(i) >= 0) {
      derivs(i) = res(i) * tau / mu;
    } else {
      derivs(i) = res(i) * (tau - 1) / mu;
    }
  }
  grad = X_t * derivs * one_over_n;
}

arma::vec huber_grad_descent_sp(const arma::colvec& y,
                                const arma::sp_mat& X,
                                const arma::sp_mat& X_t,
                                arma::vec& beta,
                                double tau, double n,
                                double one_over_n,
                                int p, int maxiter,
                                double mu,
                                double beta_tol,
                                double check_tol) {


  // gradient vector
  arma::vec grad(p);
  arma::vec last_beta = beta;
  arma::vec last_grad = grad;
  arma::vec derivs(n);

  // vector of residuals
  arma::colvec resid = y - X * beta;

  // differences from previous gradients
  // and betas
  arma::sp_vec beta_diff(p);
  arma::sp_vec grad_diff(p);

  double checkfun_diff = checkfun(resid, tau, n);
  double last_checkfun = checkfun_diff;
  double this_checkfun;

  int i = 1;
  double cross = 0;
  double delta = 1;
  // inf is the "max" norm over the vector
  while((i < maxiter) && ((arma::norm(grad, "inf") > beta_tol) || (i == 1)) &&
        (checkfun_diff * delta > check_tol || delta < 0.01)) {
    
    // step-size selection and gradient descent code originally from conquer package
    // https://github.com/XiaoouPan/conquer
    // picking step size
    delta = 1;
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(grad_diff.t() * grad_diff);
      double a2 = arma::as_scalar(beta_diff.t() * beta_diff) / cross;

      // pick the smaller and don't go _really_ off the rails
      delta = std::min(std::min(a1, a2), 1.0);
    }

    if(i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // store last gradient for step size selection
    last_grad = grad;
    update_huber_grad_sp(X_t, resid, derivs, grad, tau, mu, n, one_over_n);

    beta_diff = beta - last_beta;
    grad_diff = grad - last_grad;

    beta += (i - 1)/(i + 2) * beta_diff;
    beta += delta * grad;

    // update the beta movement for step size later
    beta_diff += delta * grad;

    // update residual vector
    resid -= X * beta_diff;
    last_beta = beta;

    this_checkfun = checkfun(resid, tau, n);
    checkfun_diff = abs((last_checkfun - this_checkfun));
    last_checkfun = this_checkfun;
    i++;
    cross = arma::as_scalar(beta_diff.t() * grad_diff);
  }
  return beta;
}


//' Compute quantile regression via accelerated gradient descent using
//' Huber approximation, warm start based on data subset
//' @param y outcome vector
//' @param X design matrix
//' @param X_sub subset of X matrix to use for "warm start" regression
//' @param y_sub subset of y to use for "warm start" regression
//' @param tau target quantile
//' @param init_beta initial guess at beta
//' @param intercept location of the intercept column, using R's indexing
//' @param num_samples number of samples used for subset of matrix used for warm start
//' @param mu neighborhood over which to smooth
//' @param maxiter maximum number of iterations to run
//' @param check_tol loss function change tolerance for early stopping
//' @param beta_tol tolerance for largest element of gradient, used
//' for early stopping
//' @param warm_start whether to start with subsamples for fitting gradient descent
//' @export
// [[Rcpp::export]]
arma::vec fit_approx_quantile_model_sp(arma::sp_mat& X,
                                       arma::colvec& y,
                                       arma::sp_mat& X_sub,
                                       arma::colvec& y_sub,
                                       double tau,
                                       arma::vec& init_beta,
                                       double mu = 1e-4,
                                       int maxiter = 10000,
                                       double beta_tol = 1e-4,
                                       double check_tol = 1e-6,
                                       int intercept = 1,
                                       double num_samples = 1000,
                                       int warm_start = 1) {

  // p is dim, n is obs, one_over_n to avoid repeated calcs
  int p = X.n_cols;
  double n = X.n_rows;

  // calc'd here because we use it a bunch
  double one_over_n = 1/n;
  // pre-transposed X
  arma::sp_mat X_t = arma::trans(X);

  // warm start
  arma::sp_mat X_t_sub = arma::trans(X_sub);

  double one_over_num_samples = 1/num_samples;

  if(warm_start == 1) {
    init_beta = huber_grad_descent_sp(y_sub,
                                      X_sub,
                                      X_t_sub,
                                      init_beta,
                                      tau, num_samples,
                                      one_over_num_samples, p,
                                      100, mu,
                                      beta_tol, check_tol);
  }

  // full data
  arma::vec beta = huber_grad_descent_sp(y, X, X_t, init_beta,
                                         tau, n, one_over_n, p,
                                         maxiter, mu,
                                         beta_tol, check_tol);


  return(beta);
}
