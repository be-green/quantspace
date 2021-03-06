% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_helpers.R
\name{getColNums}
\alias{getColNums}
\title{Get column numbers given starting values and regression specification}
\usage{
getColNums(start_list, reg_spec_data, alpha, j)
}
\arguments{
\item{start_list}{starting values (can be NA's) to be fd into sfn_start_val function}

\item{reg_spec_data}{result of ensureSpecRank function; regression matrix with full rank}

\item{alpha}{column vector of quantiles to be estimated}

\item{j}{index of quantile currently being calculated}
}
\value{
If start_list is supplied, then returns the correct column numbers
to be used in regression. Otherwise, it returns NULL.
}
\description{
Get column numbers given starting values and regression specification
}
