// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
using namespace Rcpp;


//' Drop Colinear Columns from dense matrix
//' @param X matrix to drop colinear columns from
//' @param tol tolerance for when pivot is zero in rank calculations
//' @export
// [[Rcpp::export]]
Rcpp::List qr_drop_colinear_columns(Eigen::Map<Eigen::MatrixXd>& X,
                                    const double tol = 0.000000001) {

  typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXd> CPivQR;
  typedef CPivQR::PermutationType Permutation;
  CPivQR PQR(X);

  // set tolerance for rank calculations
  // it was too low in some cases before leading
  // to inconsistent behavior
  PQR.setThreshold(tol);

  const Permutation Pmat(PQR.colsPermutation());
  const int p = X.cols();
  const int r(PQR.rank());

  if(r == X.cols()) {
    return Rcpp::List::create(
      _["spec_mat"] = X,
      _["drop_cols"] = NULL
    );
  }

  // from estimatr
  // Get all column indices
  Eigen::ArrayXi Pmat_indices = Pmat.indices();
  // Get the order for the columns you are keeping
  Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
  // Get the indices for columns you are discarding
  Eigen::ArrayXi Pmat_toss = Pmat_indices.tail(p - r);

  // copy for now so that we don't completely overwrite stuff
  Eigen::MatrixXd Xout = X;

  // move non-colinear columns up so that we can drop the end pieces
  for (Eigen::Index i=0; i < Pmat_toss.size(); i++) {
    if (Pmat_toss(i) < X.cols())
      Xout.block(0, Pmat_toss(i), X.rows(), X.cols() - Pmat_toss(i) - 1) = X.rightCols(X.cols() - Pmat_toss(i) - 1);
  }

  // just return the non-colinear columns
  return Rcpp::List::create(
    _["spec_mat"] = Xout.block(0, 0, X.rows(), Pmat_keep.size()),
    _["drop_cols"] = Pmat_toss + 1
    );
}


//' Drop Colinear Columns from sparse matrix
//' @param X matrix to drop colinear columns from
//' @param tol tolerance for when pivot is zero in rank calculations
//' @export
// [[Rcpp::export]]
Rcpp::List sparse_qr_drop_colinear_columns(Eigen::Map<Eigen::SparseMatrix<double>>& X,
                                           const double tol = 0.000000001) {
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> PQR(X);
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

  // set tolerance for rank calculations
  // it was too low in some cases before leading
  // to inconsistent behavior
  PQR.setPivotThreshold(tol);

  const int p = X.cols();
  const int r(PQR.rank());

  if(r == X.cols()) {
    return Rcpp::List::create(
      _["spec_mat"] = X,
      _["drop_cols"] = NULL
    );
  }

  // from estimatr
  // Get all column indices
  Eigen::ArrayXi Pmat_indices = Pmat.indices();
  // Get the order for the columns you are keeping
  Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
  // Get the indices for columns you are discarding
  Eigen::ArrayXi Pmat_toss = Pmat_indices.tail(p - r);

  // copy for now so that we don't completely overwrite stuff
  Eigen::MatrixXd Xout = X;

  // move non-colinear columns up so that we can drop the end pieces
  for (Eigen::Index i=0; i < Pmat_toss.size(); i++) {
    if (Pmat_toss(i) < X.cols())
      Xout.block(0, Pmat_toss(i), X.rows(), X.cols() - Pmat_toss(i) - 1) = X.rightCols(X.cols() - Pmat_toss(i) - 1);
  }

  // just return the non-colinear columns
  return Rcpp::List::create(
    _["spec_mat"] = Xout.block(0, 0, X.rows(), Pmat_keep.size()),
    _["drop_cols"] = Pmat_toss + 1
  );
}

