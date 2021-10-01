
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

#' Returns the indices of the columns to remove to construct a full rank matrix
#' @param m Rank-deficient M by N space matrix where M >= N in dgCMatrix format.
#' @param TOL tolerance for column matching
#' @returns Indices of columns to remove from m so that the remaining matrix is full rank.
#' @details Given a rank-deficient M by N sparse matrix, where M >= N, in dgCMatrix format,
#' returns the indices of columns to remove from the original matrix so that the
#' resulting matrix is full rank.
findRedundantCols = function(m, TOL = 0.000000001) {
  decomp <- qr(m)
  orig_col_names <- colnames(m)
  R <- qr.R(decomp)
  R_col_names <- colnames(R)
  if(is.null(R_col_names)) {
    R_col_names <- paste0("X", 1:ncol(R))
  }
  which(abs(diag(R)) < TOL)
}

#' Grab the last colinear column in a matrix
#' @param x matrix to drop columns from
#' @importFrom utils tail
# taken from https://stackoverflow.com/questions/12304963/using-eigenvalues-to-test-for-singularity-identifying-collinear-columns
findLastRedundantCol <- function(x) {
  xtx <- crossprod(x)

  ee <- eigen(xtx)
  evals <- zapsmall(ee$values)
  evecs <- split(zapsmall(ee$vectors),col(ee$vectors))

  cols = mapply(function(val,vec) {
    if (val!=0) NULL else which(vec!=0)
  },zapsmall(ee$values),evecs)
  l = Filter(f = function(cols) !is.null(cols), cols)
  l = unlist(utils::tail(l, 1))
  utils::tail(l, 1)
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

  init_p = ncol(spec_mat)
  p = init_p
  r = getRank(spec_mat)
  # Check if input is already matrix full rank
  if (p == r) {
    return(list(
      "spec_matrix" = spec_mat,
      "var_names" = col_names))
  }

  drop_cols = c()
  while (r < p) {
    drop_col = findLastRedundantCol(spec_mat)
    spec_mat <- spec_mat[,-drop_col]
    r = getRank(spec_mat)
    p = ncol(spec_mat)
    drop_cols = c(drop_cols, drop_col)
  }
  if(length(drop_cols) > 0) {
    nm = col_names
    if(is.null(nm)) {
      nm = 1:init_p
    }
    warning("Dropping column(s) ", paste0(nm[drop_cols], collapse = ", "),
            " due to colinearity." )
  }
  return(list(
    "spec_matrix" = spec_mat,
    "var_names" = col_names[-drop_cols]))
}

#' Convert matrix to a SparseM csr matrix
#' @param m a matrix in standard, dense format
#' @return A matrix in SparseM csr format
denseMatrixToSparse = function(m) {

  if (length(m) <= .Machine$integer.max) {
    return(as.matrix.csr(m))
  }

  # as.matrix.csr cannot coerce a long vector to csr,
  # so break the input into maximally-sized chunks with
  # are coercible and then bind them into one sparse matrix

  chunk_size <- floor(.Machine$integer.max / ncol(m)) * ncol(m)
  num_chunks <- floor(length(m) / chunk_size)

  items_to_cut <- chunk_size * num_chunks
  sparse_chunks <- lapply(
    split(t(m)[1:items_to_cut], sort(rep(1:num_chunks, chunk_size))),
    FUN = function(x) {
      return(as.matrix.csr(matrix(x, ncol=ncol(m), byrow=TRUE)))
    })

  combined <- do.call(rbind, sparse_chunks)

  if (length(m) > items_to_cut) {
    tail <- as.matrix.csr(
      matrix(
        t(m)[(items_to_cut + 1):length(m)],
        ncol=ncol(m),
        byrow=TRUE))
    combined <- rbind(combined, tail)
  }

  return(combined)
}

#' Add missing columns for specification
#' @param df data.frame data frame with data for the specification
#' @param names names of columns to include in new data.frame
#' @details Takes in a data.frame and a set of column names,
#' returns a data.frame with the specified columns, assigning 0 to the
#' value of missing columns for all rows.
addMissingSpecColumns = function(df, names) {
  missing_cols <- setdiff(names, colnames(df))
  df[missing_cols] <- 0
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

