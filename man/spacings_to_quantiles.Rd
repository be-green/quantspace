% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant_regress.R
\name{spacings_to_quantiles}
\alias{spacings_to_quantiles}
\title{Compute quantiles given parameter coefficients and data}
\usage{
spacings_to_quantiles(spacingCoef, data, jstar)
}
\arguments{
\item{spacingCoef}{J by p matrix; row is number of variables, p is number of quantiles}

\item{data}{independent variables}

\item{jstar}{index of median quantiles}
}
\value{
N by p matrix of quantiles
}
\description{
Compute quantiles given parameter coefficients and data
}
