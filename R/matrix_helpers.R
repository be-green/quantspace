
#' Function for calculating Pseudo-R^2 values
#' Evaluates the check objective function (possibly weighted) for QREG
#' @param u vector of residuals
#' @param tau probability index parameter (quantile of interest)
#' @param weight_vec optional vector of weights
#' @return Sum of absolute residuals
rho <- function(u,tau=.5,weight_vec = NULL){

  if (is.null(weight_vec)){
    ( u*(tau - (u < 0)) )
  } else {
    ( weight_vec*u*(tau - (u < 0)) )
  }

}

#' Computes the rank of a matrix. Outputs consistent results across
#' various matrix types.
#' @param m A matrix, m
#' @param TOL tolerance for rank calculation
#' @return The rank of the matrix m
getRank = function(m, TOL = 0.000000001) {

  #
  # Args:
  #   A matrix m.
  #
  # Returns:
  #   The rank of m.
  #

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
  return(setdiff(
    1:length(orig_col_names) ,
    match(
      R_col_names[which(abs(diag(R)) > TOL)],
      orig_col_names)))
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
ensureSpecFullRank = function(spec_mat, col_names) {

  # Check if input is already matrix full rank
  if (getRank(spec_mat) == ncol(spec_mat)) {
    return(list(
      "spec_matrix" = spec_mat,
      "var_names" = col_names))
  }

  transp_prod <- t(spec_mat) %*% spec_mat
  # First check if any dummy columns are all zero
  zero_cols <- which(diag(transp_prod) == 0)
  if (length(zero_cols) > 0) {
    spec_mat <- spec_mat[,-zero_cols]
    col_names <- col_names[-zero_cols]
  }

  # Check if updated matrix is full rank
  if (getRank(spec_mat) == ncol(spec_mat)) {
    return(list(
      "spec_matrix" = spec_mat,
      "var_names" = col_names))
  }

  # Use a generic routine to identify columns to drop
  dgc_spec_mat <- as(spec_mat, "dgCMatrix")
  colnames(dgc_spec_mat) <- col_names
  cols_to_drop = findRedundantCols(dgc_spec_mat)
  return(list(
    "spec_matrix" = spec_mat[,-cols_to_drop],
    "var_names" = col_names[-cols_to_drop]))
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
#' alpha: column vector of quantiles to be estimated
#' j: index of quantile currently being calculated
#' @return If start_list is supplied, then returns the correct column numbers
#' to be used in regression. Otherwise, it returns NULL.
getColNums = function(start_list,
                      reg_spec_data,
                      alpha,
                      j){

  if(!any(is.na(start_list))){
    #get the columns appropriate for starting values
    cols = reg_spec_data$var_names
    cols = paste(alpha[j], cols, sep = '_')
    col_nums = colnames(start_list) %in% cols
    return(col_nums)
  }
  return(NULL)
}
