% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxillary.R
\name{spline_weights}
\alias{spline_weights}
\alias{splw}
\title{Spline weights}
\usage{
spline_weights(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL, K, design,
  parallel = FALSE)

splw(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL, K)
}
\arguments{
\item{y}{list of observations}

\item{t}{list of time points}

\item{warp_fct}{warping function}

\item{w}{warp parameters}

\item{Sinv}{precision matricees}

\item{basis_fct}{Basis function}

\item{weights}{weights (optional)}

\item{K}{dimension}

\item{design}{design}

\item{parallel}{Calculate the sum elements in parallel?}
}
\value{
A matrix with spline weights

A matrix with spline weights
}
\description{
Spline weights

Spline weights
}
\details{
Increasing splines might possibly not work together with designs.

splw is used when no design is supplied.
splw.d is used when a desigen is supplied.
splw has been updated to handle increasing splines.
}
