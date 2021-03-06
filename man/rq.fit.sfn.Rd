% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant_regress.R
\name{rq.fit.sfn}
\alias{rq.fit.sfn}
\title{Version that complies with more general requirements}
\usage{
rq.fit.sfn(
  X,
  y,
  tau = 0.5,
  weights = NULL,
  rhs = (1 - tau) * c(t(X) \%*\% rep(1, length(y))),
  control,
  lambda,
  ...
)
}
\arguments{
\item{X}{structure of the design matrix X stored in csr format}

\item{y}{response vector}

\item{tau}{target quantile}

\item{weights}{optional vector of weights}

\item{rhs}{right hand size of dual problem}

\item{control}{control parameters for fitting routines}

\item{lambda}{ignored}

\item{...}{other arguments, ignored}
}
\description{
Version that complies with more general requirements
}
