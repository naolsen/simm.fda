### En anden ramme til krydskorrelationer hvori der kan ekstrapoleres variansmatricer.

#' dekomp
#' 
#' Positive definite square root of matrix
#'
#' @param Y 
#'
#' @return Y^{1/2}
#'
dekomp <- function(Y) {
    e <- eigen(Y)
    ev <- e$vec
    ev %*% diag(sqrt(e$val)) %*% t(ev)
}




#' Interpolate matrices
#'
#' @param t evaluation points. t[1] should correspond to a1 and t[length(t)] to a2.
#' @param a1, a2  (positive definite) matrices of same dimensions.
#'
#' @description Linearly interpolates matrices a1 and a2. Used in dynamical covariance structures.
#' Typically positive definite, but not a requirement for this function.
#' @details a1 and a2 should have same size. This is NOT checked.
#' @return An array of 
#' @export
#' @seealso kovMat
#'  
#' @examples A1 <- matrix(c(2,1,0,1,2,0.5, 0, 0.5, 2), 3, 3)
#' A2 <- diag(3)
#' A1
#' A2
#' MatInterpolate(A1, A2)
MatInterpolate <- function(t, a1, a2) {
    t <- (head(t[-1], -1) - min(t))/diff(range(t))
    A <- sapply(t, function(t) (1 - t) * a1 + t * a2)
    dim(A) <- c(dim(a1), length(t))
    A
}

# MatInterpolate(ti, mat, mat + 1)[, , ]


#' Dynamical covariance structure
#'
#' @param t evalutation points
#' @param a list of matrices
#' @param tw Anchor points. Must be in increasing order
#' @param timefunction Temporal covariance structure
#' @param stack How should stack the observations. TRUE corresponds to the implementation of ppMUlti, while FALSE is the order ordering.
#' @param noise Noise. Either a single number or a vector of length K*length(t)
#' @param ... Parameters passed to timefunction
#' 
#' @details 
#'  
#'
#' @return A positive definite covariance matrix
#' @export
#'
#' @examples
#' @references See pdf simm.fda for introduction and of details of dynamical covariance structure
kovMat <- function(t, a, tw, timefunction, stack = TRUE, noise = 0, ...) {
    z <- length(tw)
    adim <- ncol(a[[1]])
    lt <- length(t)
    
    bmat <- matrix(nc = adim * lt, nr = adim)
    # tw2 <- tw[-1]
    for (j in 1:(z - 1)) {
        g <- t >= tw[j] & t <= tw[j + 1]
        ti <- t[g]
        if (all(!g)) next
        h <- MatInterpolate(c(tw[j], ti, tw[j + 1]), a[[j]], a[[j + 1]])

        bmat[, rep(g, each = adim)] <- (apply(MatInterpolate(c(tw[j], ti, tw[j + 1]), a[[j]], a[[j + 1]]), 
            3, dekomp))
    }
    
    if (stack) {
        bmat <- bmat[, as.vector(t(matrix(1:(adim * length(t)), nr = adim)))]
        return((t(bmat) %*% bmat) * (matrix(1, nc = adim, nr = adim) %x% outer(t, t, timefunction, ...)) + 
            diag(x = noise, lt * adim))
    }
    
    (t(bmat) %*% bmat) * (outer(t, t, timefunction, ...) %x% matrix(1, nc = adim, nr = adim)) + diag(x = noise, 
        lt * adim)
    # outer(bmat, bmat, FUN='%*%')
    
}

# dim(kovMat(ti, list(mat, mat + 1), c(0, 1), function(x, y) { abs(x - y) }, stack = TRUE))




#' Temporal covariance structures
#'
#' @param s,t  time values 
#' @param lambda drift
#' @rdname Tidsfunktioner
#' 
#' @export
#'
#' @examples
OUtid <- function(s, t, lambda) { ## OU process
    exp(-abs(s - t) * lambda)
}

## tidsfunktion for Brownsk Bro & Bev?gelse.

#'
#'
#' @param bridge Brownian bridge or Brownian motion?
#' @details Time covariance functions should not be called directly. Either use as part of kovMat or call using outer.
#'
#' @rdname Tidsfunktioner
#'
#' @export
#'
#' @examples
BMtid <- function(s, t, bridge = TRUE) { ## Brown bridge/motion
    x <- pmin(s, t)
    if (bridge) return(x - s*t)
    else return(x)
}


## Matern temporal covariance structure


#'
#' @param range range parameter
#' @param smooth smoothness parameter
#' @description Matern covariance structure: 
#' 2^(smooth-1)/gamma(smooth) * range^smooth ^ K_smooth(|s-t|/range)
#' 
#' @rdname Tidsfunktioner
#' @return
#' @export
#'
#' @examples
Materntid <- function(s,t, range, smooth) {

  range <-abs(s-t)/range
  range[range == 0] <- 1e-12
  (besselK(range, smooth) * range^smooth / (gamma(smooth) *2 ^(smooth-1)))
}



## tidsfunktion for Matern kovarians med bro.
Matern.bro <- function(s,t, range, smooth) {
  x <- s*t
  range <-abs(s-t)/range
  range[range == 0] <- 1e-12
  (besselK(range, smooth) * range^smooth / (gamma(smooth) *2 ^(smooth-1))) * (pmin(s, t) - x)
}



#' Polynomial bridge
#'
#' @param tmin, tmax Max and min values
#' @param coef Coefficients, i.e. a marginal variance given by coef[1] + t*(1-t)*coef[2]
#'
#' @details Polynomial bridge cannot be used as a timefunction. It is recommended to combine this with timefunction.bridge
#' or something like it to ensure the right structure.
#' 
#'
#' @return If done correctly, a 'bridge structure covariance', see
#' @export
#'
#' @seealso timefunction.bridge
poly.bridge <- function(tmin, tmax, coef) { ## Bemærk ændringer på parametriseringen ift. tidligere kode.
  coef[2]*tmin*(1-tmax) + coef[1]
}

#' Polynomial 'bridge' of arbitrary order
#'
#' @param tmin 
#' @param tmax 
#' @param coef Polynomial coefficients for bridge, in increasing order. 
#' 
#' @details Generalise \code{poly.bridge} to arbitrary order.
#'
#' @return
#' @export
#'
#' @seealso \link{poly.bro}, timefunction.bridge
poly.larger <- function(tmin, tmax, coef) {
  
  tminmax <- tmin*(1-tmax)
  mat <- coef[2]*tminmax + coef[1]
  if (length(coef) > 2) for (i in 3:length(coef))
  {
    mat <- mat + coef[i] * tminmax^i
  }
  mat
}


#' Matern covariance and bridge
#'
#' @param s,t 
#' @param t 
#' @param range 
#' @param smooth 
#' @param koef 
#' 
#' @details This is a specific implementation for Matern covariance.
#'
#' @return
#' @export
#'
#' @seealso \link{timefunction.bridge} 
#' 
#' @examples 
#' f <- function(s, t, params) poly.Matern(s, t, params[1], params[2], params[3:4] )
#' 
#' a1 <- diag(c(1.5, 0.1, 0.7))
#' a2 <- a1 + 2
#' 
#' ti <- seq(0,1, 0.02)
#' 
#' kovMat(ti, list(a1, a2), c(0,1), f, params = c(0.5, 2, 1, 1), noise = 0.1)
#' 
poly.Matern <- function(s,t, range, smooth, koef) {
  x1 <- pmin(s,t)
  x2 <- pmax(s,t)
  range <- (x2-x1)/range
  range[range == 0] <- 1e-12
  (besselK(range, smooth) * range^smooth / (gamma(smooth) *2 ^(smooth-1))) * poly.bridge(x1, x2, koef)
}


#' Timefunction and bridge
#'
#' Creates a 'bridged' time function, i.e. a combination of an 'ordinary' time function and a bridge covariance. 
#' As the result is a time function, it can be used in conjunction with kovMat.
#'
#' @param s,t 
#' @param coef Coefficients for bridge
#' @param timefunction A timefunction
#' @param ... Parameters passed to timefunction
#'
#' @return A time function structure.
#' @export
#'
#' @examples 
#' 
#' f <- function(s, t, param) timefunction.bridge(s, t, coef = c(param[1],param[2]), OUtid, lambda = param[3])
#' 
#' ti <- seq(0, 1, 0.02)
#' 
#' outer(ti , ti, f, param = c(0.2, 1, 0.7)) 
#' outer(ti , ti, f, param = c(1, 0, 0.7)) 
#' # Try make a contour plot to see the differences
#' 
#' 
timefunction.bridge <- function(s, t, coef, timefunction,  ...) {
  x1 <- pmin(s,t)
  x2 <- pmax(s,t)
  timefunction(s, t, ...) * poly.bridge(x1, x2, coef)
}


