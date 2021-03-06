% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Misc.R
\name{make_dyn_cov}
\alias{make_dyn_cov}
\title{Parametrize dynamical covariance structure using Cholesky decomposition}
\usage{
make_dyn_cov(timefunction, knots, K)
}
\arguments{
\item{timefunction}{Time function for dynamical covariance}

\item{knots}{knots for dynamical covariance}

\item{K}{Dimension}
}
\value{
A function
}
\description{
Cholesky decomposition is the recommended way of parametrising the positive definite matrices used for dynamical covariance structure.
This function provides a handy wrapper for doing this.
}
\details{
.
}
\examples{

g <- make_dyn_cov(Materntid, c(0, 0.5, 1), 3)

# assume we have some parameters in a vector pars, e.g. pars = c(rep(c(1,0,0, 1,0,1), 3), 0.1, 1.5) (diagonal covariance)
g(ti, param = pars[1:18], range = pars[19], smooth = pars[20], noise = 1)

}
\seealso{
\link{kovMat}
}
