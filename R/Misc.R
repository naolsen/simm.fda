




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

#' Bounded optimization using nlm
#'
#' @param fct Function to be optimized over. One vector argument.
#' @param p Starting values
#' @param lower lower bounds for paramters.
#' @param upper upper bounds for parameters.
#' @param ... Further arguments to nlm
#' 
#' @description Performs bounded (box-constrained) non-linear minimization using the fast but unconstrained \code{nlm} optimizer. 
#' This function is in many cases been faster than box-constrained optimization using \link{optim}.
#'
#' @details 
#' nlm.bound doesn't check validity of the bounds, but does check if starting values satisfy the constraints.
#' If any starting value is not strictly between its prescribed bounds, it is adjusted with a warning. 
#'
#' @return
#'
#' @seealso \link{nlm}
nlm.bound <- function(fct, p, lower, upper, ...) {
  
  pl <- length(p)
  if (length(upper) != pl || length(lower) != pl) stop("lengths of p, lower and upper must match!")
  
  for (i in 1:pl) {
    if (p[i] <= lower[i]) {
      p[i] <- 0.999* lower[i]+ 0.001* upper[i]
      warning("p[i] is on or outside of the lower boundary. Adjusting by 0.001* (upper-lower)")
    }
    else if (p[i] >= upper[i]) {
      p[i] <- 0.001* lower[i]+ 0.999* upper[i]
      warning("p[i] is on or outside of the upper boundary. Adjusting by -0.001* (upper-lower)")
    }
  }
  
  nul <- p
  
  nedre <- p - lower
  ovre <- upper - p
  p0 <- rep(0, length(p))
  
  p <-
    nlm(f = function(par) {

      p2 <- exp(-abs(par))
      par <- nul + (ovre* (par > 0) - nedre* (par < 0))*(1-p2)
      
      fct(par)
    } , p = p0, ...)
  p2 <- exp(-abs(pest <- p$estimate))
  
  #pr <- sign(p$estimate)
  
  p$estimate <- nul + (ovre* (pest > 0) - nedre* (pest < 0))*(1-p2)
  
  p
}


