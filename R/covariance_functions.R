

#make_warp_cov <- function(covariance_fct, default.values) {
#  attr(covariance_fct, "param") <- c(default.values)
#  covariance_fct
#}

#' Brownian Bridge warp covariance
#'
#' @param default.value Default value
#' 
#' @description This is called default covariance function as this is typically used
#' @return
#' @export
#'
#' @examples
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
#' @return
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
#' @return
#' @export
#'
shift_covariance <- function(default.value) {
  g <- function(t, param) {
    param^2
  }
  attr(g, 'param') <- c(shift = default.value)
  g
}
  
  


