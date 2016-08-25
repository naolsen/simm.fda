




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
#' @description Performs bounded non-linear minimization using the fast but unconstrained nlm optimizer.
#'
#' @details nlm.bound doesn't check validity of bounds, but does check if starting values satisfy the constraints.
#' If any starting value is not strictly between its prescribed bounds, it is adjusted with a warning. 
#' 
#'
#' @return output from nlm, but with estimates transformed back.
#' 
#'
#' @seealso \link{nlm}
nlm.bound <- function(fct, p, lower, upper, ...) {
  
  pl <- length(p)
  if (length(upper) != pl || length(lower) != pl) stop("lengths of p, lower and upper must match!")
  
  for (i in pl) {
    if (p[i] <= lower[i]) {
      p[i] <- 0.999* lower+ 0.001* upper
      warning("p[i] is on or outside of the lower boundary. Adjusting by 0.001* (upper-lower)")
    }
    else if (p[i] >= upper[i]) {
      p[i] <- 0.001* lower+ 0.999* upper
      warning("p[i] is on or outside of the upper bound. Adjusting by -0.001* (upper-lower)")
    }
  }
  p0 <- (p - lower) / (upper - lower)
  p0 <- log(p0 / (1- p0))
  
  p <-
    nlm(f = function(par) {
      par <- exp(par)
      par[is.infinite(par)] <- 1e300
      par <- lower + par / (1+par) *(upper - lower)
      #par <- lower + exp(par)/(1+exp(par))*(upper-lower)
      fct(par)
    } , p = p0, ...)
  p$estimate <- lower + exp(p$estimate)/ ( 1+ exp(p$estimate)) *(upper - lower)
  p
}

