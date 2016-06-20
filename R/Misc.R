




#' Parametrize dynamical covariance structure using Cholesky decomposition
#'
#' @param timefunction Time function for dynamical covariance
#' @param knots knots for dynamical covariance
#' @param K Dimension 
#' 
#' 
#' @description Cholesky decomposition is the recommended way of parametrising the positive definite matrices used for dynamical covariance structure.
#' This function provides a handy wrapper for doing this.
#'
#' @details .
#' @return A function
#' @export
#' @seealso \link{kovMat}
#'
#' @examples
#' 
#' g <- make_dyn_cov(Materntid, c(0, 0.5, 1), 3)
#' 
#' # assume we have some parameters in a vector pars, e.g. pars = c(rep(c(1,0,0, 1,0,1), 3), 0.1, 1.5) (diagonal covariance)
#' g(ti, param = pars[1:18], range = pars[19], smooth = pars[20], noise = 1)
#' 
make_dyn_cov <- function( timefunction, knots, K) {
  
  indices <- which(lower.tri(matrix(0,K,K), diag = T))
  lpar <- length(indices)
  d <- length(knots)
  cat("This function will have", lpar*d, "parameters \n")
  
  return( function(t, param, ...) {
    a <- list()
    for (i in 1:d) {
      mat <- matrix(0, K, K)
      mat[indices] <- param[((d-1)*lpar+1):(d*lpar) ]
      a[[i]] <- mat %*% t(mat)
    }
    knots
    kovMat(t, a, knots, timefunction, ...)
    
    
  })
}
