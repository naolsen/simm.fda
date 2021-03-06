
#' Jacobian of warped B-spline function
#'
#' Computes the Jacobian of a warped B-spline function.
#' @param t (warped) evaluation points.
#' @param dwarp Jacobian of the warping function for the given evaluation points.
#' @param basis_fct basis function.
#' @param c spline weights.
#' 
#' 

Zi <- function(t, dwarp, basis_fct, c) {
  basis <- basis_fct(t, deriv = TRUE)
  dwarp <- dwarp * ((basis %*% c)[ , 1])
  return(dwarp)
}


## Multivariat version af Zi

multi.Zi <- function(t, dwarp, basis_fct, c, mw) {
  K <- ncol(c)
  b0 <- basis_fct(t, deriv = TRUE)
  
  basis <- array( , dim = c(length(t), K, mw))
  for (i in 1:K) basis[,i,] <- dwarp* as.numeric(b0 %*% c[,i])
  dim(basis) <- c(length(t)*K, mw)
  
  return(basis)
}
