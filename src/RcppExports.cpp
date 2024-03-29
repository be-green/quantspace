// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// qr_drop_colinear_columns
Rcpp::List qr_drop_colinear_columns(Eigen::Map<Eigen::MatrixXd>& X, const double tol);
RcppExport SEXP _quantspace_qr_drop_colinear_columns(SEXP XSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(qr_drop_colinear_columns(X, tol));
    return rcpp_result_gen;
END_RCPP
}
// sparse_qr_drop_colinear_columns
Rcpp::List sparse_qr_drop_colinear_columns(Eigen::Map<Eigen::SparseMatrix<double>>& X, const double tol);
RcppExport SEXP _quantspace_sparse_qr_drop_colinear_columns(SEXP XSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double>>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_qr_drop_colinear_columns(X, tol));
    return rcpp_result_gen;
END_RCPP
}
// checkfun
double checkfun(arma::vec& res, double tau);
RcppExport SEXP _quantspace_checkfun(SEXP resSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type res(resSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(checkfun(res, tau));
    return rcpp_result_gen;
END_RCPP
}
// parallelVectorCheckFun
double parallelVectorCheckFun(arma::vec& x, double tau);
RcppExport SEXP _quantspace_parallelVectorCheckFun(SEXP xSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelVectorCheckFun(x, tau));
    return rcpp_result_gen;
END_RCPP
}
// reorder_columns
void reorder_columns(arma::mat& X, int intercept);
RcppExport SEXP _quantspace_reorder_columns(SEXP XSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    reorder_columns(X, intercept);
    return R_NilValue;
END_RCPP
}
// fast_rexp
arma::vec fast_rexp(int n);
RcppExport SEXP _quantspace_fast_rexp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_rexp(n));
    return rcpp_result_gen;
END_RCPP
}
// update_huber_grad
void update_huber_grad(const arma::mat& X_t, const arma::vec& res, arma::vec& derivs, arma::vec& grad, double tau, double mu, int n, double one_over_n);
RcppExport SEXP _quantspace_update_huber_grad(SEXP X_tSEXP, SEXP resSEXP, SEXP derivsSEXP, SEXP gradSEXP, SEXP tauSEXP, SEXP muSEXP, SEXP nSEXP, SEXP one_over_nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X_t(X_tSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type one_over_n(one_over_nSEXP);
    update_huber_grad(X_t, res, derivs, grad, tau, mu, n, one_over_n);
    return R_NilValue;
END_RCPP
}
// z_score
void z_score(arma::mat& X, const arma::rowvec& colwise_avg_x, const arma::vec& colwise_sd_x, const int p);
RcppExport SEXP _quantspace_z_score(SEXP XSEXP, SEXP colwise_avg_xSEXP, SEXP colwise_sd_xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type colwise_avg_x(colwise_avg_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type colwise_sd_x(colwise_sd_xSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    z_score(X, colwise_avg_x, colwise_sd_x, p);
    return R_NilValue;
END_RCPP
}
// arma_qr_drop_colinear_columns
Rcpp::List arma_qr_drop_colinear_columns(arma::mat& X);
RcppExport SEXP _quantspace_arma_qr_drop_colinear_columns(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_qr_drop_colinear_columns(X));
    return rcpp_result_gen;
END_RCPP
}
// huber_grad_descent
arma::vec huber_grad_descent(const arma::colvec& y, const arma::mat& X, const arma::mat& X_t, arma::vec& beta, arma::vec& grad, arma::vec& derivs, double tau, double n, double one_over_n, int p, int maxiter, double mu, double beta_tol, double check_tol, double min_delta);
RcppExport SEXP _quantspace_huber_grad_descent(SEXP ySEXP, SEXP XSEXP, SEXP X_tSEXP, SEXP betaSEXP, SEXP gradSEXP, SEXP derivsSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP one_over_nSEXP, SEXP pSEXP, SEXP maxiterSEXP, SEXP muSEXP, SEXP beta_tolSEXP, SEXP check_tolSEXP, SEXP min_deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_t(X_tSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type one_over_n(one_over_nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tol(beta_tolSEXP);
    Rcpp::traits::input_parameter< double >::type check_tol(check_tolSEXP);
    Rcpp::traits::input_parameter< double >::type min_delta(min_deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(huber_grad_descent(y, X, X_t, beta, grad, derivs, tau, n, one_over_n, p, maxiter, mu, beta_tol, check_tol, min_delta));
    return rcpp_result_gen;
END_RCPP
}
// fit_approx_quantile_model
Rcpp::List fit_approx_quantile_model(arma::mat& X, arma::vec& y, arma::mat& X_sub, arma::vec& y_sub, double tau, arma::colvec init_beta, double mu, int maxiter, double beta_tol, double check_tol, int intercept, double num_samples, int warm_start, int scale, double lambda, double min_delta);
RcppExport SEXP _quantspace_fit_approx_quantile_model(SEXP XSEXP, SEXP ySEXP, SEXP X_subSEXP, SEXP y_subSEXP, SEXP tauSEXP, SEXP init_betaSEXP, SEXP muSEXP, SEXP maxiterSEXP, SEXP beta_tolSEXP, SEXP check_tolSEXP, SEXP interceptSEXP, SEXP num_samplesSEXP, SEXP warm_startSEXP, SEXP scaleSEXP, SEXP lambdaSEXP, SEXP min_deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_sub(X_subSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y_sub(y_subSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type init_beta(init_betaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tol(beta_tolSEXP);
    Rcpp::traits::input_parameter< double >::type check_tol(check_tolSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< double >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< int >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type min_delta(min_deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_approx_quantile_model(X, y, X_sub, y_sub, tau, init_beta, mu, maxiter, beta_tol, check_tol, intercept, num_samples, warm_start, scale, lambda, min_delta));
    return rcpp_result_gen;
END_RCPP
}
// fit_penalize_approx_quantile_model
arma::vec fit_penalize_approx_quantile_model(arma::mat& X, arma::vec& y, arma::mat& X_sub, arma::vec& y_sub, double tau, arma::colvec init_beta, double mu, int maxiter, double beta_tol, double check_tol, int intercept, double num_samples, int warm_start, int scale);
RcppExport SEXP _quantspace_fit_penalize_approx_quantile_model(SEXP XSEXP, SEXP ySEXP, SEXP X_subSEXP, SEXP y_subSEXP, SEXP tauSEXP, SEXP init_betaSEXP, SEXP muSEXP, SEXP maxiterSEXP, SEXP beta_tolSEXP, SEXP check_tolSEXP, SEXP interceptSEXP, SEXP num_samplesSEXP, SEXP warm_startSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_sub(X_subSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y_sub(y_subSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type init_beta(init_betaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tol(beta_tolSEXP);
    Rcpp::traits::input_parameter< double >::type check_tol(check_tolSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< double >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< int >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_penalize_approx_quantile_model(X, y, X_sub, y_sub, tau, init_beta, mu, maxiter, beta_tol, check_tol, intercept, num_samples, warm_start, scale));
    return rcpp_result_gen;
END_RCPP
}
// glob_obs_mat
arma::mat glob_obs_mat(const arma::mat& X, const arma::vec& r, double thresh);
RcppExport SEXP _quantspace_glob_obs_mat(SEXP XSEXP, SEXP rSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(glob_obs_mat(X, r, thresh));
    return rcpp_result_gen;
END_RCPP
}
// glob_obs_vec
arma::vec glob_obs_vec(const arma::vec& y, const arma::vec& r, double thresh);
RcppExport SEXP _quantspace_glob_obs_vec(SEXP ySEXP, SEXP rSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(glob_obs_vec(y, r, thresh));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_quantspace_qr_drop_colinear_columns", (DL_FUNC) &_quantspace_qr_drop_colinear_columns, 2},
    {"_quantspace_sparse_qr_drop_colinear_columns", (DL_FUNC) &_quantspace_sparse_qr_drop_colinear_columns, 2},
    {"_quantspace_checkfun", (DL_FUNC) &_quantspace_checkfun, 2},
    {"_quantspace_parallelVectorCheckFun", (DL_FUNC) &_quantspace_parallelVectorCheckFun, 2},
    {"_quantspace_reorder_columns", (DL_FUNC) &_quantspace_reorder_columns, 2},
    {"_quantspace_fast_rexp", (DL_FUNC) &_quantspace_fast_rexp, 1},
    {"_quantspace_update_huber_grad", (DL_FUNC) &_quantspace_update_huber_grad, 8},
    {"_quantspace_z_score", (DL_FUNC) &_quantspace_z_score, 4},
    {"_quantspace_arma_qr_drop_colinear_columns", (DL_FUNC) &_quantspace_arma_qr_drop_colinear_columns, 1},
    {"_quantspace_huber_grad_descent", (DL_FUNC) &_quantspace_huber_grad_descent, 15},
    {"_quantspace_fit_approx_quantile_model", (DL_FUNC) &_quantspace_fit_approx_quantile_model, 16},
    {"_quantspace_fit_penalize_approx_quantile_model", (DL_FUNC) &_quantspace_fit_penalize_approx_quantile_model, 14},
    {"_quantspace_glob_obs_mat", (DL_FUNC) &_quantspace_glob_obs_mat, 3},
    {"_quantspace_glob_obs_vec", (DL_FUNC) &_quantspace_glob_obs_vec, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_quantspace(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
