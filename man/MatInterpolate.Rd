% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kovMat.R
\name{MatInterpolate}
\alias{MatInterpolate}
\title{Interpolate matrices}
\usage{
MatInterpolate(t, a1, a2)
}
\arguments{
\item{t}{evaluation points. t[1] should correspond to a1 and t[length(t)] to a2.}

\item{a1, }{a2  (positive definite) matrices of same dimensions.}
}
\value{
A three-dimensional array of interpolated matrices
}
\description{
Linearly interpolates matrices a1 and a2. Used in dynamical covariance structures.
Typically positive definite, but not a requirement for this function.
}
\details{
a1 and a2 should have same size. This is NOT checked.
}
\examples{
A1 <- matrix(c(2,1,0,1,2,0.5, 0, 0.5, 2), 3, 3)
A2 <- diag(3)
A1
A2
MatInterpolate(A1, A2)
}
\seealso{
\link{kovMat}
}
