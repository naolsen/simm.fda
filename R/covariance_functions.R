

#make_warp_cov <- function(covariance_fct, default.values) {
#  attr(covariance_fct, "param") <- c(default.values)
#  covariance_fct
#}

#' Brownian Bridge warp covariance
#'
#' @param default.value Default value
#' 
#' @description This is called default covariance function as this is typically used
#' @return Covariance function
#' @export
#'
#' @examples warp_cov <- default_warp_cov(0.5)
#' warp_cov
#' 
default_warp_cov <- function(default.value) {
  g <- function(t, param) {
    Brown(t, c(param^2, 0 ), motion = FALSE)
  }
  attr(g, 'param') <- c(tau = default.value)
  g
}

#' Brownian Bridge warp covariance + shift parameter
#'
#' @param default.values c(shift variance, Brownian Bridge variance)
#'
#' @description Combines default_warp_cov and a shift parameter. Shift and warp are assumed to be independent.
#' 
#' @return Covariance function
#' @export
#'
#' @seealso \link{default_warp_cov}, \link{w.shift}
#' @examples
warp_and_shift_cov <- function(default.values) {
  g <- function(t, param) {
    
    m <- length(t)
    wc <- matrix(0,nc= m+1, nr= m+1)
    wc[1] <- param[1]^2
    wc[2:(m+1), 2:(m+1)] <- Brown(t, c(param[2]^2, 0 ), motion = FALSE)

    wc
  }
  attr(g, 'param') <- c(shift = default.values[1], tau = default.values[2])
  g
  
}

#' Simple one-parameter covariance for pure shift
#'
#' @param default.value Default value
#'
#' @return Covariance function
#' @export
#'
shift_covariance <- function(default.value) {
  g <- function(t, param) {
    param^2
  }
  attr(g, 'param') <- c(shift = default.value)
  g
}
  
#' Second-order derivative covariance
#'
#' @param default.value Default value
#' @param type If bridge, use integrated brownian bridge (one parameter), else use integrated brownian motion (two parameters)
#'
#' @return Covariance function
#' @export
#'
order2_covariance <- function(default.value, type ="bridge") {

   if (type == "bridge") g <- function(t, param)  
      param^2* outer(t,t, function(x,y) x*y*(1/3 + 1/6*(x^2+y^2)-1/2*pmax.int(x,y)) - 1/6* pmin.int(x,y)^3)

   else g <- function(t, param)
      outer(t,t, function(x,y) param[1]^2*x*y +  param[2]^2* pmin.int(x,y)^2*(pmax.int(x,y) - pmin.int(x,y) /3))
   
  attr(g, 'param') <- c(tau = default.value)
  #attr(g, 'inv') <- invers.order2

  g
}



## Experimental
if (F) {

invers.order2 <- local({
  result <- NA
  t.last <- -1
  
  f <- function(t, param) {
    if (!all(t == t.last)) 
     {
      res <- outer(t,t, function(x,y) x*y*(1/3 + 1/6*(x^2+y^2)-1/2*pmax.int(x,y)) - 1/6* pmin.int(x,y)^3)
      invers <- chol2inv(chol(res))
      invers[abs(invers) < max(invers)* 1e-5] <- 0
      result <<- as(invers, "dgCMatrix")
      t.last <<- t
    }
    (1/param^2)*result
  }
})
}


