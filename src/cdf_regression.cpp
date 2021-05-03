// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::plugins(cpp11)]]

//' Calculate full cdf regression via accelerated gradient descent using
//' Huber approximation
//' @param y outcome vector
//' @param X design matrix
//' @param min_tau smallest tau to fit
//' @param max_tau largest tau to fit
//' @param tau_step step size between min and max tau
//' @param init_beta initial guess at beta
//' @param mu neighborhood over which to smooth
//' @param maxiter maximum number of iterations to run
//' @param check_tol loss function change tolerance for early stopping
//' @param beta_tol tolerance for largest element of gradient, used
//' for early stopping
arma::mat cdf_regression(arma::mat& X,
                         arma::colvec& y,
                         double min_tau,
                         double max_tau,
                         double tau_step,
                         arma::colvec& init_beta,
                         double mu = 0.00001,
                         int maxiter = 1000,
                         double check_tol = 1e-4,
                         double beta_tol = 1e-4) {

  // dim
  int p = X.n_cols;

  // obs, double for division
  double n = X.n_rows;
  double one_over_n = 1/n;
  double init_tau = 0.5;
  double delta = 2;

  int num_taus = floor((max_tau - init_tau)/tau_step +
                       (init_tau - min_tau)/tau_step + 1);
  int median_col = floor((init_tau - min_tau)/tau_step);

  arma::mat beta_mat(p, num_taus);

  arma::colvec resid(n);
  arma::colvec init_coef = huber_grad_descent(X, y, init_tau, init_beta, mu,
                                              maxiter, check_tol, beta_tol);
  arma::colvec coef_guess = init_coef;
  arma::mat X_t = arma::trans(X);
  arma::colvec pointwise_derivs(n);
  // arma::colvec tau_grad(p);
  arma::colvec median_resid = y - X * init_coef;
  arma::colvec real_coef(p);

  Rcpp::Rcout << "beta_mat has : " << beta_mat.n_cols << " columns.\n";
  Rcpp::Rcout << "beta_mat has : " << beta_mat.n_rows << " rows.\n";

  double new_tau = init_tau + tau_step;
  resid = median_resid;
  int i = median_col;
  beta_mat.col(i) = init_coef;
  arma::vec init_tau_grad = dbeta_dtau(y, X, real_coef, init_tau,
                                       mu, n, one_over_n);
  arma::vec tau_grad = init_tau_grad;
  // increase tau relative to 0.5
  while(i < (num_taus - 1)) {

    Rcpp::checkUserInterrupt();

    i = i + 1;
    // approximate new coefficients by moving along the gradient
    // of the loss function
    Rcpp::Rcout<< "Tau is " << new_tau << "\n";
    Rcpp::Rcout<< "Tau gradient is " << sum(tau_step * tau_grad) << "\n";

    // get good guess of optimal coefs w/ first order approx of loss
    // function by taking derivative w/r/t tau
    // then get real optimal coefficients via accelerated gradient descent
    coef_guess = coef_guess + tau_step * tau_grad;

    real_coef = huber_grad_descent(X, y, new_tau, coef_guess, mu,
                                   maxiter, check_tol, beta_tol);
    resid = y - X * real_coef;

    // update guess for new stuff
    coef_guess = real_coef;

    // fill relevant column with updated coefficients
    beta_mat.col(i) = real_coef;
    tau_grad = dbeta_dtau(y, X, real_coef, new_tau,
                          mu, n, one_over_n);
    new_tau = new_tau + tau_step;
  }

  i = median_col;
  resid = median_resid;
  coef_guess = init_coef;
  tau_step = -tau_step;
  new_tau = init_tau + tau_step;
  tau_grad = init_tau_grad;
  // increase tau relative to 0.5
  while(i > 0) {
    i = i - 1;
    // // approximate new coefficients by moving along the gradient
    // // of the loss function
    // // we are decreasing tau now

    Rcpp::Rcout<< "Tau gradient is " << sum(tau_step * tau_grad) << "\n";

    // // get good guess of optimal coefs w/ first order approx of loss
    // // function by taking derivative w/r/t tau
    // // then get real optimal coefficients via accelerated gradient descent
    coef_guess = coef_guess + tau_step * tau_grad;

    Rcpp::checkUserInterrupt();

    real_coef = huber_grad_descent(X, y, new_tau, coef_guess, mu,
                                   maxiter, check_tol, beta_tol);
    coef_guess = real_coef;
    resid = y - X * real_coef;
    // fill relevant column with updated coefficients
    beta_mat.col(i) = real_coef;
    tau_grad = dbeta_dtau(y, X, real_coef, new_tau,
                          mu, n, one_over_n);
    new_tau = new_tau + tau_step;
  }
  return beta_mat;
}


// gradient w/r/t quantile tau
double dh_dtau(double res, double mu) {
  if (res > mu) {
    return(res - mu/2);
  } else if(res < -mu) {
    return(-res - mu/2);
  } else {
    return((res * res) / (2 * mu));
  }
}

arma::colvec dh_dtau(arma::colvec x, double mu) {
  arma::colvec v(x.n_elem);
  for(int i = 0; i < x.n_elem; i = i + 1) {
    v(i) = dh_dtau(x(i), mu);
  }
  return(v);
}

arma::colvec dbeta_dtau(arma::colvec& y, arma::mat& X, arma::colvec& beta,
                        double tau, double mu, int n, double one_over_n) {

  int p = X.n_cols;
  arma::colvec resid = y - X * beta;
  double dl_dtau = mean(dh_dtau(resid, mu));
  Rcpp::Rcout<< "Loss function gradient is " << dl_dtau << "\n";
  Rcpp::Rcout<< "Tau is " << tau << "\n";

  double sum_resid_gt_0 = 0;
  double sum_resid_lt_0 = 0;

  arma::colvec x_gt_0(p);
  arma::colvec x_lt_0(p);

  x_gt_0.zeros(p);
  x_lt_0.zeros(p);
  double num_pos_resids = 0;
  double num_neg_resids = 0;

  for(int i=0; i < n; i = i + 1) {
    if(resid(i) >= 0) {
      num_pos_resids += 1;
      sum_resid_gt_0 += resid(i);
      x_gt_0 = x_gt_0 + X.row(i).t();

    } else {

      num_neg_resids += 1;
      sum_resid_lt_0 += resid(i);
      x_lt_0 = x_lt_0 + X.row(i).t();
    }
  }

  double pos_normalizer = 1/num_pos_resids;
  double neg_normalizer = 1/num_neg_resids;

  arma::colvec first_denom = neg_normalizer * (tau - 1) * x_lt_0;
  arma::colvec second_denom = pos_normalizer * tau * x_gt_0;

  Rcpp::Rcout << "Avg Pos Resid: " << pos_normalizer * sum_resid_gt_0 << " \n";
  Rcpp::Rcout << "Avg Neg Resid: " << neg_normalizer * sum_resid_lt_0 << " \n";
  Rcpp::Rcout << "Bottom Colvec: " << first_denom + second_denom << " \n";
  Rcpp::Rcout << "Inv Colvec: " << 1/(first_denom + second_denom)  << " \n";



  return (dl_dtau - pos_normalizer * sum_resid_gt_0 -
          neg_normalizer * sum_resid_lt_0) /
            (first_denom + second_denom);
}
