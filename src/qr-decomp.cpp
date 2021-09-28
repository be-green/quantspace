// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat get_Q (arma::mat& X) {
  arma::mat Q;
  arma::mat R;

  arma::qr_econ(Q, R, X);

  arma::umat P_mat;

  qr(Q, R, P_mat, X, "matrix");
  return Q;
}

// [[Rcpp::export]]
int get_rank (arma::mat& X) {
    return arma::rank(X);
}
