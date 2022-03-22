// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List qr_drop_colinear_columns(Eigen::Map<Eigen::MatrixXd>& X) {

  typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXd> CPivQR;
  typedef CPivQR::PermutationType Permutation;
  const CPivQR PQR(X);
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

// [[Rcpp::export]]
Rcpp::List sparse_qr_drop_colinear_columns(Eigen::Map<Eigen::SparseMatrix<double>>& X) {
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> PQR(X);
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

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

