






double.integrated.warp <- function(tw) {
  mw <- length(tw) + 1
  tw <- c(0, tw, 1)
  if (any(diff(tw) <= 0)) stop("tw must be an increasing sequence between 0 and 1.")
  dtw <- diff(tw)
  dtw2 <- diff(tw^2) / 2
  
  gradmat <- matrix(0, mw+1, mw)
  for (i in 1:mw) { ## derivatives for 
    w <- numeric(mw)
    w[i] <- 1
    vt <- tw* c(0,cumsum(dtw*w)) - c(0,cumsum(dtw2*w))
    
    gradmat[,i] <- vt- tw*vt[length(vt)]
  }
  gradmat <- gradmat[2:mw, ] ## initial and ending zeros not needed
  
  v <- function(w, t, w_grad = FALSE) {
    
    if (!w_grad) {
      stopifnot(length(w) == mw)
      #vt <- tw* c(0,cumsum(dtw*w)) - c(0,cumsum(dtw2*w))
      #t +approx(tw, vt- tw*vt[length(vt)], xout = t, rule = 2)$y
      
      vt <- tw* c(0,cumsum(dtw*w)) - c(0,cumsum(dtw2*w))
      vt <- tw + vt- tw*vt[length(vt)]
      if (any(diff(vt) <= 0)) vt <- tw
      approx(tw, vt, xout = t, rule = 2)$y
    }
    else {
      # Derivative of warp function
      # Note: does not depend on w
      apply(cbind(tw[1:(mw-1)], tw[2:mw], tw[3:(mw+1)]), 1, function(x) {
        a <- rep(0, length(t))
        a[t > x[1] & t < x[2]] <- ((t - x[1])/(x[2] - x[1]))[t > x[1] & t < x[2]]
        a[t >= x[2] & t < x[3]] <- (1 - ((t - x[2])/(x[3] - x[2])))[t >= x[2] & t < x[3]]
        return(a)
      }) %*% gradmat
    }
  }
  attr(v, 'initialize') <- 0
  attr(v, 'mw') <- mw
  
  return(v)
  
}



