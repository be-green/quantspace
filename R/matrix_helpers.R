
#' Function for calculating Pseudo-R^2 values
#' Evaluates the check objective function (possibly weighted) for QREG
#' @param u vector of residuals
#' @param tau probability index parameter (quantile of interest)
#' @param weights optional vector of weights
#' @return Sum of absolute residuals
rho <- function(u,tau=.5,weights = NULL){

  if (is.null(weights)){
    ( u*(tau - (u < 0)) )
  } else {
    ( weights*u*(tau - (u < 0)) )
  }

}

#' Computes the rank of a matrix. Outputs consistent results across
#' various matrix types.
#' @param m A matrix, m
#' @param TOL tolerance for rank calculation
#' @return The rank of the matrix m
getRank = function(m, TOL = 1e-10) {
  transp_prod <- as.matrix(t(m) %*% (m))
  return(sum(abs(diag(qr.R(qr(transp_prod)))) > TOL))
}


#' Ensure that a regression specification is full rank
#' @details Verifies if a regression specification is full-rank. If the input
#' is rank-deficient, identifies and drops columns so that the remaining
#' matrix is full-rank.
#' @param spec_mat An M by N regression specification matrix, where M > N.
#' @param col_names Column names of the regression specification.
#' @return A full-rank specification matrix. If the input is full-rank, returns the
#' input unmodified. Otherwise, returns a matrix with a subset of the columns
#' from the input.
#' @importFrom methods as
ensureSpecFullRank = function(spec_mat, col_names) {
  init_names <- col_names
  p <- length(init_names)
  if(is.matrix(spec_mat)) {
    spec_mat <- qr_drop_colinear_columns(spec_mat)
  } else {
    spec_mat <- csrToDgc(spec_mat)
    spec_mat <- sparse_qr_drop_colinear_columns(spec_mat)
    spec_mat$spec_mat <- matrixTocsr(spec_mat$spec_mat)
  }

  if(length(spec_mat$drop_cols) > 1 || spec_mat$drop_cols > 0) {
    col_names <- col_names[-spec_mat$drop_cols]
  }

  r <- length(col_names)

  if(r < p) {
    warning("Dropping ", p - r, " colinear columns: ",
             paste0(init_names[spec_mat$drop_cols], collapse = ", "), ".")
  }

  return(list(
    "spec_matrix" = spec_mat$spec_mat,
    "var_names" = col_names))
}

#' @import methods
#' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

#' Convert Matrix Sparse Matrix to SparseM row compressed matrix
#' @param X column compressed matrix from Matrix package
#' @importFrom SparseM as.matrix.csr
#' @importFrom methods new
#' @importClassesFrom SparseM matrix.csc matrix.csr
matrixTocsc = function(X) {
  X_csc <- methods::new("matrix.csc", ra = X@x,
               ja = X@i + 1L,
               ia = X@p + 1L,
               dimension = X@Dim)
  X_csc
}

#' Convert Matrix Sparse Matrix to SparseM row compressed matrix
#' @param X column compressed matrix from Matrix package
#' @importFrom SparseM as.matrix.csr
#' @importFrom methods new
#' @importClassesFrom SparseM matrix.csc matrix.csr
matrixTocsr = function(X) {
  X_csc <- methods::new("matrix.csc", ra = X@x,
               ja = X@i + 1L,
               ia = X@p + 1L,
               dimension = X@Dim)
  SparseM::as.matrix.csr(X_csc)
}

#' Convert SparseM row compressed matrix to Matrix dgC matrix
#' @param X column compressed matrix from Matrix package
#' @importFrom SparseM as.matrix.csr
#' @importFrom Matrix sparseMatrix
csrToDgc = function(X) {
  Matrix::sparseMatrix(
    p = X@ia - 1L,
    j = X@ja,
    x = X@ra
  )
}

#' Convert matrix to a SparseM csr matrix
#' @param m a matrix in standard, dense format
#' @return A matrix in SparseM csr format
#' @importFrom Matrix Matrix
denseMatrixToSparse = function(m) {
  m = Matrix::Matrix(m, sparse = T)
  matrixTocsr(m)
}

#' Add missing columns for specification
#' @param df data.frame data frame with data for the specification
#' @param names names of columns to include in new data.frame
#' @details Takes in a data.frame and a set of column names,
#' returns a data.frame with the specified columns, assigning 0 to the
#' value of missing columns for all rows.
addMissingSpecColumns = function(df, names) {
  missing_cols <- setdiff(names, colnames(df))
  df[missing_cols] <- NA
  return(df[names])
}

#' Get column numbers given starting values and regression specification
#' @param start_list starting values (can be NA's) to be fd into sfn_start_val function
#' @param reg_spec_data result of ensureSpecRank function; regression matrix with full rank
#' @param alpha column vector of quantiles to be estimated
#' @param j index of quantile currently being calculated
#' @return If start_list is supplied, then returns the correct column numbers
#' to be used in regression. Otherwise, it returns NULL.
getColNums = function(start_list,
                      reg_spec_data,
                      alpha,
                      j){

  if(!any(is.na(start_list))){
    cols = reg_spec_data$var_names
    cols = paste(alpha[j], cols, sep = '_')
    col_nums = colnames(start_list) %in% cols
    return(col_nums)
  }
  return(NULL)
}

#' For copying matrices as in Matlab (works for sparse matrices)
#' @param X matrix to replicate
#' @param m number of times replicate the matrix rows
#' @param n number of times replicate the matrix columns
repMat <- function(X, m, n){
  Y <- do.call(rbind, rep(list(X), m))
  do.call(cbind, rep(list(Y), n))
}

#' Create sparse diagonal matrix with vector x on diagonal
#' @importFrom SparseM diag
#' @importFrom SparseM as.matrix.csr
#' @param v Vector to use as the diagonal of the matrix
spDiag <- function(v){
  SparseM::as.matrix.csr(SparseM::diag(v))
}

#' Return column sums of matrix
#' @param m matrix to sum up
spSums <- function(m){
  N <- dim(m)[1]
  ones <- denseMatrixToSparse(repMat(1,1,N))
  return(denseMatrixToSparse(ones%*%m))
}

#' Inverse of a matrix, but catches the error
#' @param a matrix to invert
inv <- function (a)
{
  if (length(a) == 0)
    return(matrix(0, nrow = 0, ncol = 0))
  if ((!is.numeric(a) && !is.complex(a)) || !is.matrix(a))
    stop("Argument 'a' must be a numeric or complex matrix.")
  if (nrow(a) != ncol(a))
    stop("Matrix 'a' must be square.")
  e <- try(b <- solve(a), silent = TRUE)
  if (inherits(e, "try-error")) {
    warning("Matrix appears to be singular.")
    b <- rep(Inf, length(a))
    dim(b) <- dim(a)
  }
  return(b)
}

