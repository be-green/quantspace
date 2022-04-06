// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
# include <RcppParallel.h>
# include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::plugins(cpp11)]]

double checkfun (double x, double tau) {
  return x * (tau - (x < 0));
}

// [[Rcpp::export]]
double checkfun (arma::vec& res, double tau) {
  double value = 0;
  int n = res.n_elem;
  for(int i = 0; i < n; i++) {
    value -= checkfun(res(i), tau);
  }
  return value;
}


// This is taken from the RcppParallel tutorial
// https://rcppcore.github.io/RcppParallel/#parallelreduce
// Edited to work with an arma vec + check function
struct SumCheckFun : public RcppParallel::Worker
{
  // source vector
  const arma::vec& input;

  // target quantile
  double tau;

  // accumulated value
  double value;


  // constructors
  SumCheckFun(const arma::vec& input, double tau) : input(input), tau(tau), value(0) {}
  SumCheckFun(const SumCheckFun& sum, RcppParallel::Split) : input(sum.input), tau(sum.tau), value(0) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    value += std::accumulate(input.begin() + begin,
                             input.begin() + end, 0.0,
                             [&](double x, double y) {return x + checkfun(y, tau);});
  }

  // join avg_y value with that of another SumCheckFun
  void join(const SumCheckFun& rhs) {
    value += rhs.value;
  }
};

// [[Rcpp::export]]
double parallelVectorCheckFun(arma::vec& x, double tau) {

  // declare the SumCheckFun instance
  SumCheckFun scf(x, tau);

  // call parallel_reduce to start the work
  RcppParallel::parallelReduce(0, x.n_elem, scf);

  // return the computed sum
  return scf.value;
}


// [[Rcpp::export]]
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
void z_score(arma::mat& X, const arma::rowvec& colwise_avg_x, const arma::vec& colwise_sd_x, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - colwise_avg_x(i)) / colwise_sd_x(i);
  }
}

// [[Rcpp::export]]
Rcpp::List arma_qr_drop_colinear_columns(arma::mat& X) {
  arma::mat Q;
  arma::mat R;
  arma::uvec P;

  int p = X.n_cols;

  // qr decomposition
  qr(Q, R, P, X, "vector");

  int r = arma::rank(R);

  Rcpp::Rcout << "P rows: "<< P.n_rows << "\n";
  Rcpp::Rcout << "P cols: "<< P.n_cols << "\n";

  if(r == p) {
    Rcpp::List lout = Rcpp::List::create(
      Rcpp::Named("spec_mat") =  X,
      Rcpp::Named("drop_cols") = 0
    );
    return(lout);
  }

  arma::uvec drop_cols = P.rows(r - 1, p - 1);
  arma::uvec keep_cols = P.rows(0, r - 2);

  X = X.cols(keep_cols);
  Rcpp::List lout = Rcpp::List::create(
    Rcpp::Named("spec_mat") =  X,
    Rcpp::Named("drop_cols") = drop_cols + 1
    );

  return(lout);

}

// [[Rcpp::export]]
arma::vec huber_grad_descent(const arma::colvec& y, const arma::mat& X,
                             const arma::mat& X_t, arma::vec& beta,
                             arma::vec& grad, arma::vec& derivs,
                             double tau, double n, double one_over_n,
                             int p, int maxiter, double mu,
                             double beta_tol, double check_tol,
                             double min_delta = 1e-10) {

  // gradient vector
  arma::vec last_beta = beta;
  arma::vec last_grad = grad;

  // vector of residuals
  arma::vec resid = y - X * beta;

  // differences from previous gradients
  // and betas
  arma::vec beta_diff(p);
  beta_diff.zeros();
  arma::vec grad_diff(p);

  // hyperparams for backtracking line search
  // idk how to optimize these but for now they are pretty
  // arbitrary
  double c = 0.5;
  double reduction_factor = 0.8;
  double loss = parallelVectorCheckFun(resid, tau);
  double t = sum(c * grad);
  double loss_diff = -100;
  double last_loss = loss;

  arma::vec resid_update(n);

  double i = 1;
  double init_delta = 1 / std::max(tau, 1 - tau);
  double delta = init_delta;
  int exitflag = 0;

  // inf is the "max" norm over the vector
  while((i < maxiter) && ((arma::norm(grad, "inf") > beta_tol) || (i == 1)) &&
        abs(loss_diff / delta) > check_tol && exitflag == 0) {


    if(std::fmod(i, 100) == 0) {
      Rcpp::checkUserInterrupt();
    }

    update_huber_grad(X_t, resid, derivs, grad, tau, mu, n, one_over_n);

    // acceleration

    // 0 in the first step
    beta += (i - 1)/(i + 2) * beta_diff;

    // gradient update from projection
    beta += delta * grad;

    resid = y - X * beta;

    loss = parallelVectorCheckFun(resid, tau);

    loss_diff = last_loss - loss;
    t = sum(c * grad);

    // if we don't reduce loss function enough, make step size smaller
    while(loss_diff < delta * t && exitflag == 0)  {
      Rcpp::checkUserInterrupt();

      if(delta < min_delta) {
        exitflag = 1;
      }

      beta -= delta * grad;
      delta *= reduction_factor;

      // move projection a single gradient step
      beta += delta * grad;

      resid = y - X * beta;
      loss = parallelVectorCheckFun(resid, tau);
      loss_diff = last_loss - loss;
    }
    // if it's all good, update the beta
    last_beta = beta;
    last_loss = loss;
    i++;
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
//' @param lambda optional lasso penalty weight
//' @param min_delta smallest allowed step size for gradient descent
//' @export
// [[Rcpp::export]]
Rcpp::List fit_approx_quantile_model(arma::mat& X,
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
                                     int scale = 1,
                                     double lambda = 0,
                                     double min_delta = 1e-10) {

  // p is dim, n is obs, one_over_n to avoid repeated calcs
  int p = X.n_cols;
  double n = X.n_rows;

  arma::vec q = {tau};
  // calc'd here because we use it a bunch
  double one_over_n = 1/n;

  arma::rowvec colwise_avg_x;
  arma::vec colwise_sd_x;
  double avg_y;
  double sd_y;

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
    colwise_avg_x = arma::mean(X, 0);
    colwise_sd_x = arma::stddev(X, 0, 0).t();

    // standardize X
    for (int i = 0; i < X.n_cols; i++) {
      X.col(i) = (X.col(i) - colwise_avg_x(i)) / colwise_sd_x(i);

      if(warm_start == 1) {
        X_sub.col(i) = (X_sub.col(i) - colwise_avg_x(i)) / colwise_sd_x(i);
      }
    }

    avg_y = arma::mean(y);
    sd_y = arma::stddev(y);
    if(intercept > 0) {
      y -= avg_y;
      if(warm_start == 1) {
        y_sub -= avg_y;
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

  if (lambda > 1e-8) {
    arma::mat R(X.n_cols, X.n_cols, arma::fill::zeros);
    R.diag() += lambda;

    R.rows(arma::span(1, X.n_cols - 1));
    arma::vec r(p);
    r.zeros();

    // add positive component of lasso penalty
    X = arma::join_vert(X, R * n);
    y = arma::join_vert(y, r);

    // add negative component of lasso penalty
    X = arma::join_vert(X, R * -n);
    y = arma::join_vert(y, r);

    X_sub = arma::join_vert(X_sub, R * num_samples);
    y_sub = arma::join_vert(y_sub, r);

    X_sub = arma::join_vert(X_sub, R * -num_samples);
    y_sub = arma::join_vert(y_sub, r);

    num_samples += r.n_rows * 2;

    n = n + R.n_rows * 2;
  }

  // pre-transposed X
  arma::mat X_t = arma::trans(X);

  // initialize w. 0 gradient
  arma::vec grad(p);
  grad.zeros();

  if(warm_start == 1) {

    // warm start
    arma::mat X_t_sub = arma::trans(X_sub);

    double one_over_num_samples = 1/num_samples;

    arma::vec derivs(num_samples);
    derivs.zeros();

    init_beta = huber_grad_descent(y_sub,
                                   X_sub,
                                   X_t_sub,
                                   init_beta,
                                   grad, derivs,
                                   tau, num_samples,
                                   one_over_num_samples, p,
                                   5, mu,
                                   beta_tol, check_tol, min_delta);

  }


  init_beta(0) = arma::as_scalar(arma::quantile(y - X.cols(1, p - 1) * init_beta.rows(1, p - 1),
                                 q));

  arma::vec residuals = y - X * init_beta;
  arma::vec derivs(n);

  update_huber_grad(X_t, residuals, derivs, grad, tau, mu, n, one_over_n);

  // full data
  arma::vec beta =  huber_grad_descent(y, X, X_t, init_beta,
                                       grad, derivs,
                                       tau, n, one_over_n, p,
                                       maxiter, mu,
                                       beta_tol, check_tol,min_delta);

  if(lambda > 1e-8) {
    n = n - (2 * p);
    X = X.rows(arma::span(0, n - 1));
    y = y.rows(arma::span(0, n - 1));
  }

  if(scale == 1) {

    // unstandardize
    y += avg_y;
    double m;
    double s;

    // unstandardize coefficients
    if(intercept > 0) {

      for (int i = 1; i < X.n_cols; i++) {
        m = colwise_avg_x(i - 1);
        s = colwise_sd_x(i - 1);
        X.col(i) = (X.col(i) + m) * s;
      }
      beta.rows(1, p - 1) /= colwise_sd_x;
      beta(0) += avg_y - arma::as_scalar(colwise_avg_x * beta.rows(1, p - 1));

    } else {

      for (int i = 0; i < p; i++) {
        m = colwise_avg_x(i);
        s = colwise_sd_x(i);
        X.col(i) = (X.col(i) + m) * s;
      }

      beta.rows(0, p - 1) /= colwise_sd_x;
    }
    // if the intercept column of X isn't the first column,
    // re-order the coefficients
    if(intercept > 1) {
      reorder_columns(X, intercept);
      beta.insert_rows(intercept - 1, beta(0));
      beta.shed_row(0);
    }
  }

  arma::vec fitted_vals = X * beta;
  return(Rcpp::List::create(
      Rcpp::Named("coefficients") = beta,
      Rcpp::Named("residuals") = y - fitted_vals,
      Rcpp::Named("fitted") = fitted_vals
  ));
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
arma::vec fit_penalize_approx_quantile_model(arma::mat& X,
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

  arma::rowvec colwise_avg_x;
  arma::vec colwise_sd_x;
  double avg_y;
  double sd_y;

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
    colwise_avg_x = arma::mean(X, 0);
    colwise_sd_x = arma::stddev(X, 0, 0).t();

    // standardize X
    for (int i = 0; i < X.n_cols; i++) {
      X.col(i) = (X.col(i) - colwise_avg_x(i)) / colwise_sd_x(i);

      if(warm_start == 1) {
        X_sub.col(i) = (X_sub.col(i) - colwise_avg_x(i)) / colwise_sd_x(i);
      }
    }

    avg_y = arma::mean(y);
    sd_y = arma::stddev(y);
    if(intercept > 0) {
      y -= avg_y;
      if(warm_start == 1) {
        y_sub -= avg_y;
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


  // initialize w. 0 gradient
  arma::vec grad(p);
  grad.zeros();

  if(warm_start == 1) {

    // warm start
    arma::mat X_t_sub = arma::trans(X_sub);

    double one_over_num_samples = 1/num_samples;

    arma::vec derivs(num_samples);
    derivs.zeros();


    init_beta = huber_grad_descent(y_sub,
                                   X_sub,
                                   X_t_sub,
                                   init_beta,
                                   grad, derivs,
                                   tau, num_samples,
                                   one_over_num_samples, p,
                                   100, mu,
                                   beta_tol, check_tol);

  }



  init_beta(0) = arma::as_scalar(arma::quantile(y - X.cols(1, p - 1) * init_beta.rows(1, p - 1),
                                 q));

  arma::vec residuals = y - X * init_beta;
  arma::vec derivs(n);

  update_huber_grad(X_t, residuals, derivs, grad, tau, mu, n, one_over_n);

  // full data
  arma::vec beta =  huber_grad_descent(y, X, X_t, init_beta,
                                       grad, derivs,
                                       tau, n, one_over_n, p,
                                       maxiter, mu,
                                       beta_tol, check_tol);
  if(scale == 1) {

    // unstandardize
    y += avg_y;
    double m;
    double s;

    // unstandardize coefficients
    if(intercept > 0) {
      for (int i = 1; i < X.n_cols; i++) {
        m = colwise_avg_x(i - 1);
        s = colwise_sd_x(i - 1);
        X.col(i) = (X.col(i) + m) * s;
      }
      beta.rows(1, p - 1) /= colwise_sd_x;
      beta(0) += avg_y - arma::as_scalar(colwise_avg_x * beta.rows(1, p - 1));
    } else {

      for (int i = 0; i < p; i++) {
        m = colwise_avg_x(i);
        s = colwise_sd_x(i);
        X.col(i) = (X.col(i) + m) * s;
      }

      beta.rows(0, p - 1) /= colwise_sd_x;
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


//' Glob observations w/ residuals above a certain magnitude
//' @param X design matrix
//' @param r vector of residuals
//' @param thresh threshold to use when computing globs
//' @export
// [[Rcpp::export]]
arma::mat glob_obs_mat(const arma::mat& X, const arma::vec& r, double thresh) {
  arma::uvec keep_indices = arma::find(r < thresh && r > -thresh);
  arma::uvec glob_indices_upper = arma::find(r > thresh);
  arma::uvec glob_indices_lower = arma::find(r < -thresh);
  arma::mat globbed_x_upper = sum(X.rows(glob_indices_upper), 0);
  arma::mat globbed_x_lower = sum(X.rows(glob_indices_lower), 0);

  arma::mat x_globbed = arma::join_vert(
    globbed_x_upper,
    X.rows(keep_indices),
    globbed_x_lower
  );
  return x_globbed;
}

//' Glob observations w/ residuals above a certain magnitude
//' @param y design matrix
//' @param r vector of residuals
//' @param thresh threshold to use when computing globs
// [[Rcpp::export]]
arma::vec glob_obs_vec(const arma::vec& y, const arma::vec& r, double thresh) {

  arma::uvec keep_indices = arma::find(r < thresh && r > -thresh);
  arma::uvec glob_indices_upper = arma::find(r > thresh);
  arma::uvec glob_indices_lower = arma::find(r < -thresh);
  arma::vec globbed_y_upper = sum(y.rows(glob_indices_upper));
  arma::vec globbed_y_lower = sum(y.rows(glob_indices_lower));

  arma::mat y_globbed = arma::join_vert(
    globbed_y_upper,
    y.rows(keep_indices),
    globbed_y_lower
  );
  return y_globbed;
}
