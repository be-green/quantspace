% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_helpers.R
\name{findRedundantCols}
\alias{findRedundantCols}
\title{Returns the indices of the columns to remove to construct a full rank matrix}
\usage{
findRedundantCols(m, TOL = 1e-09)
}
\arguments{
\item{m}{Rank-deficient M by N space matrix where M >= N in dgCMatrix format.}

\item{TOL}{tolerance for column matching}
}
\value{
Indices of columns to remove from m so that the remaining matrix is full rank.
}
\description{
Returns the indices of the columns to remove to construct a full rank matrix
}
\details{
Given a rank-deficient M by N sparse matrix, where M >= N, in dgCMatrix format,
returns the indices of columns to remove from the original matrix so that the
resulting matrix is full rank.
}
