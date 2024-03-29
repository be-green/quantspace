% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_effects.R
\name{marginal_effects}
\alias{marginal_effects}
\title{Get all marginal effects of variables in the fit}
\usage{
marginal_effects(fit, type = "mea", variable = "all", data = NA, trim = 0.05)
}
\arguments{
\item{fit}{model fitted by \code{qs()}}

\item{type}{one of "ame" (average marginal effects) or
"mea" (marginal effects at the average)}

\item{variable}{which variable to calculate marginal effects on}

\item{data}{optional data.frame that specifies level of data to calculate
marginal effects}

\item{trim}{What to trim the variable of interest at, 0 < trim < 0.5}
}
\description{
Get all marginal effects of variables in the fit
}
\details{
The trim defaults to using the 95th percentile
instead of the max because there may be large outliers. You can over-ride
by setting trim to 0, which will use the min and max.
By default, marginal effects will calculate marginal effects for
all variables.
}
