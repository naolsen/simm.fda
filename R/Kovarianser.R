### Se ogs? R-pakke 'MitRbib'

## The values are stacked such that { 1, ..., n } ## refers to the first response observed at t_1, ...,
## t_n, and { n+1, ..., 2n } refers to the second response observed at t_1, ..., t_n etc.



## New function for minimum . Actually equal to pmin
minFkt <- function(s, t) {
    x <- s > t
    x * t + (s * (!x))
}




#' Brownian motion and Brownian bridge
#'
#' @param t 
#' @param par Two-dimensional. Variance of BB/BM and variance of noise.
#' @param motion logical. Brownian bridge or brownian motion
#'
#' @return
#' @export
#'
#' @examples
Brown <- function(t, par = c(tau = 1, noise = 1), motion = FALSE) {
    S <- outer(t, t, minFkt)
    if (!motion) 
        S <- S - outer(t, t)
    
    S * par[1] + diag(x = par[2], nrow = length(t))
}

## Matern
#' Matern covariance
#' 
#' Matern covariance with noise
#'
#' @param t time points
#' @param par parameters
#'
#' @return Covariance matrix for Matern covariance
#' @export
#'
#' @examples 
#' t <- seq(0,1 , 0.01)
#' Maternk(t, c(0.1, 2, 1))
Maternk <- function( t, par = c(range =1, smooth=2, noise = 1))  {
  ran <- par[1]
  glat <- par[2]
  l <- length(t)
  S <- outer(t,t, FUN = function(x,y) { abs(x-y)})/ran
  S[S == 0] <- 1e-12
  (besselK(S, glat) * S^glat / (gamma(glat) *2 ^(glat-1))) + diag(l)*par[3]
  
}

## Multivariat blandet Matern-kovarians med st?j 1
poly.Matern.kov <- function( t, sig, range =1, smooth=200, koef2) {
  ran <- range
  glat <- smooth
  l <- length(t)
  S <- outer(t,t, FUN = function(x,y) { abs(x-y)})/ran
  S[S == 0] <- 1e-12
  #(besselK(S, glat) * S^glat / (gamma(glat) *2 ^(glat-1))) *  (1+koef2*(outer(t,t, minFkt) - outer(t, t))) + diag(l)*par[3]
  sig %x% ((besselK(S, glat) * S^glat / (gamma(glat) *2 ^(glat-1))) *  
         (1+koef2*(outer(t,t, minFkt) - outer(t, t)))) +  diag(x = 1,   l*nrow(sig))
  
}


## OU proces
#' Ornstein-Uhlenbeck process
#'
#' @param t time points
#' @param par vector with drift and noise terms.
#'
#' @details Note that the marginal variance (expect for noise) is assumed to be 1.
#' @return Covariance matrix for OU process
#' @export
#'
#' @examples
OUproces <- function(t, par = c(lambda = 1, noise = 0)) {
    S <- outer(t, t, FUN = function(x, y) {
        abs(x - y)
    })
    exp(-par[1] * S) + diag(par[2], length(t))
}




#' Multivariate OU process
#'
#' @param t 
#' @param lambda drift parameter.
#' @param sig Marginal covariance. K x K matrix 
#' @param noise 
#'
#' @return
#' @export
#'
#' @examples
mvOUproces <- function(t, lambda = 1, sig, noise = 0) {
    S <- outer(t, t, FUN = function(x, y) {
        abs(x - y)
    })
    sig %x% (exp(-lambda * S)) + diag(x = rep(noise, each = length(t)))
}

mvBrown <- function(t, tau = 1, sig, noise = 0, motion = FALSE) {
    sig %x% Brown(t, par = c(tau, 0), motion) + diag(x = rep(noise, each = length(t)))
}


#' Multivariate Matern covariance
#'
#' @param t 
#' @param range,smooth range and smoothness parameters 
#' @param smooth 
#' @param noise. noise for different coordinates. Must be a vector of length K. 
#' @param sig Marginal covariance. K x K matrix 
#' 
#' @details noise and sig should fit in dimension
#'
#' @return
#' @export
#'
#' @examples
#' 
mvMatern <- function( t, range =1, smooth=2, noise = 1, sig)  {
  sig %x% Maternk(t, c(range, smooth, 0)) +  diag(x = rep(noise, each = length(t)))
}


## Wrapper functions for optimization algorithms

## wrapper for multivariat Brownsk bro. Bruger nedre trekantsmatrix i Cholesky dekompositionen
## Skalaparameteren par[1] er redundant.
BBwrapper <- function(t, par) {
    pa <- matrix(c(par[2], par[5], par[6], 0, par[3], par[7], 0, 0, par[4]), nrow = 3)
    mvBrown(t, par[1], pa %*% t(pa), c(1, 1, 1), motion = FALSE)
}

## Multivariat OU proces
mvOUwrapper <- function(t, par) {
    pa <- matrix(c(par[2], par[5], par[6], 0, par[3], par[7], 0, 0, par[4]), nrow = 3)
    mvOUproces(t, par[1], pa %*% t(pa), c(1, 1, 1))
} 




### Smart cholesky dekomposition for ukorrelerede multi-dim processer.

c.mv.BM <- function(t,  par, motion = TRUE) {
  s <- length(t)
  z <- list()
  
  S <- outer(t, t, minFkt)
  if (!motion)     S <- S - outer(t, t)
  for (i in 1:length(par)) z[[i]] <- chol(S*par[i] + diag(s))
  bdiag(z)
  
}


c.mv.Matern <- function( t, range =1, smooth=200, par ) {
  s <- length(t)
  z <- list()
  
  S <- outer(t,t, FUN = function(x,y) { abs(x-y)})/range
  S[S == 0] <- 1e-12
  
  
  for (i in 1:length(par)) z[[i]] <-  
    chol(par[i]* (besselK(S, smooth) * S^smooth / (gamma(smooth) *2 ^(smooth-1))) + diag(s))
  
  bdiag(z)
  
}


#' Diagonal covariance
#'
#' @param t 
#' @param par parameters. One parameter for each dimension
#' 
#' @description This is the simplest covariance function. 
#'
#' @return
#' @export
#'
#' @examples 
#' # Three-dimensional curve:
#' amp_cov <- diag_covariance
#' amp_cov_par <- c(1,1,1)
#' ppMulti(y,t, etc...)
diag_covariance <- function(t, par) {
  Diagonal(x = rep(par, each = length(t)))
}

attr(diag_covariance, "inv_cov_fct") <- function(t, par) {
  Diagonal(x = rep(1/par, each = length(t)))
}

  