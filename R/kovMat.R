### En anden ramme til krydskorrelationer hvori der kan ekstrapoleres variansmatricer.

#' dekomp
#'
#' @param Y 
#'
#' @return 
#'
dekomp <- function(Y) {
    e <- eigen(Y)
    ev <- e$vec
    ev %*% diag(sqrt(e$val)) %*% t(ev)
}


# Lars: faster replacement 2x faster for 200 x 200 matrices

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
#' @param s, t time values 
#' @param lambda drift
#' @rdname Tidsfunktioner
#' 
#' @return
#' @export
#'
#' @examples
OUtid <- function(s, t, lambda) { ## OU process
    exp(-abs(s - t) * lambda)
}

## tidsfunktion for Brownsk Bro & Bev?gelse.


#' Temporal covariance structures
#'
#'
#' @param bridge Brownian bridge or Brownian motion?
#' @details Time covariance functions should not be called directly.
#'
#' @rdname Tidsfunktioner
#' @return
#' @export
#'
#' @examples
BMtid <- function(s, t, bridge = TRUE) { ## Brown bridge/motion
    x <- pmin(s, t)
    if (bridge) return(x - s*t)
    else return(x)
}


## tidsfunktion for Matern kovarians.
#' Matern temporal covariance structure
#'

#' @param range range parameter
#' @param smooth smoothness parameter
#' @description Matern covariance structure: 
#' 2^(smooth-1)/gamma(smooth) * range^smooth ^ K_smooth(|s-t|/range)
#' @details Time covariance functions should not be called directly. 
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

# 
#  (kovMat(ti, list(mat, mat + 1), c(0, 1), OUtid, stack = TRUE, lambda = 100)) 
# kovMat(ti, list(mat, mat + 1, mest), c(0, 0.5, 1), BBtid, stack = TRUE, lambda = 1, noise = 1) 
# 
# (kovMat(ti, list(mat, mat, mat +1), c(0, 0.5, 1), Matern.bro, stack = TRUE, range = 1, smooth=3)) -0
# 
# mvMatern(ti, 1, 3, c(0,0), mat)


## tidsfunktion for Matern kovarians med bro.
Matern.bro <- function(s,t, range, smooth) {
  x <- s*t
  range <-abs(s-t)/range
  range[range == 0] <- 1e-12
  (besselK(range, smooth) * range^smooth / (gamma(smooth) *2 ^(smooth-1))) * (pmin(s, t) - x)
}

## Endnu mere dynamisk seriel korrelation. Boer kombineres med en af de oevrige modeller
## a*t^2*(1-t)^2 + b*t*(1-t) + c
## Bemaerk at i praksis vil der vaere overparametrisering, saa saet en af parametrerne til en.
## Parametrisering: tmin og tmax: hhv. s \wedge t og s \vee t (ie. pmin og pmax)
poly.bro <- function(tmin, tmax, koef) {
#  x1 <- pmin(s,t)
#  x2 <- 1-pmax(s,t)
  koef[1]*tmin^2*(1-tmax)^2 + koef[2]*tmin*(1-tmax) + koef[3]
}

poly.larger <- function(tmin, tmax, koef) {
  #  x1 <- pmin(s,t)
  #  x2 <- 1-pmax(s,t)
  koef[1]*tmin^2*(1-tmax)^2 + koef[2]*tmin*(1-tmax) + koef[3]
}


poly.Matern <- function(s,t, range, smooth, koef) {
  x1 <- pmin(s,t)
  x2 <- pmax(s,t)
  range <- (x2-x1)/range
  range[range == 0] <- 1e-12
  (besselK(range, smooth) * range^smooth / (gamma(smooth) *2 ^(smooth-1))) * poly.bro(x1, x2, koef)
}

#outer(ti, ti, poly.Matern, 1, 3, par = c(1,1,1))
