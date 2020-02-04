###
### Likelihood functions


#' Linearised likelihood functions.
#'
#' @param param Amplitude parameters
#' @param param.w Warp parameters
#' @param r,Zis Residual vectors/matrices
#' @param amp_cov Amplitude covariance function
#' @param warp_cov Warp covariance function
#' @param t Observation time points
#' @param tw 
#' @param sig Return sigma or likelihood?
#' @param pr Print option. Only used for debugging.
#' @rdname likelis
#' 
#' @details \code{likelihood} is the default option. \code{like.par} performs parallelized likelihood. 
#' \code{like.nowarp} is used when there is no warping used (and thus no linearisation).
#'
#' @return -2*logL(parameters)
#' @export
#'
#' @seealso \link{ppMulti}
likelihood <- function(param, r, amp_cov, t, param.w, Zis, warp_cov, tw, sig=FALSE, pr = FALSE, w = NULL) {
  
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, param.w)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }
  
  n <- length(r)
  m <- sapply(r, length)
  
  sq <- logdet <- 0
  for (i in 1:n) {
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], param)
      
      U <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(U)) return(1e10) ## Exception handling;
    } else {
      
      S <- U <- diag(1, m[i])
    }
    rr <- as.numeric(r[[i]])
    ZZ <- Zis[[i]]
    
    if (!is.null(warp_cov)) {
      A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
      LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
      x <- t(A) %*% rr
    } else {
      LR <- x <- 0
    }
    sq <- sq + (sum(backsolve(U, rr, transpose = TRUE)^2)
                - t(x) %*% LR %*% x)
    logdet_tmp <- 0
    if (!is.null(warp_cov)) logdet_tmp <- determinant(LR)$modulus[1]
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }
  
  if (!is.null(warp_cov)) logdet <- logdet - n * determinant(Cinv)$modulus[1]
  
  sigmahat <- as.numeric(sq /sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  if (pr)  print(res)
  if (sig) return (sqrt(sigmahat))
  return(min(res, 1e10))
  
}

## Parallelized likelihood


#'
#' @rdname likelis
#' 
like.par <- function(param, r, amp_cov, t, param.w, Zis, warp_cov, tw, sig=FALSE, pr = FALSE, w = NULL) {
  
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, param.w)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }
  
  n <- length(r)
  m <- sapply(r, length)
  
  sqLogd <- foreach(i = 1:n, ZZ = Zis, rr = r, tid = t, .combine = '+', .noexport = c("Zis", "r", "t", "param_w"), .inorder = FALSE) %dopar% {
    if (!is.null(amp_cov)) {
      S <- amp_cov(tid, param)
      
      U <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(U)) return(1e10) ## Exception handling;
    } else {
      S <- U <- diag(1, m[i])
    }
    
    
    rr <- as.numeric(rr)
    
    if (!is.null(warp_cov)) {
      A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
      LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
      x <- t(A) %*% rr
    } else {
      LR <- x <- 0
    }
    sqq <- (sum(backsolve(U, rr, transpose = TRUE)^2) - t(x) %*% LR %*% x)
    
    logdet_tmp <- 0
    if (!is.null(warp_cov)) logdet_tmp <- determinant(LR)$modulus[1]
    logd <- - (logdet_tmp - 2 * sum(log(diag(U))))
    c(sqq, logd)
  }
  
  if (!is.null(warp_cov)) sqLogd[2] <- sqLogd[2] - n * determinant(Cinv)$modulus[1]
  
  sigmahat <- as.numeric(sqLogd[1] /sum(m))
  res <- sum(m) * log(sigmahat) + sqLogd[2]
  if (pr)  print(res)
  if (sig) return (sqrt(sigmahat))
  return(min(res, 1e10))
  
}   


## Used when no warp is given
#'
#' @rdname likelis
#' 
like.nowarp <- function(param,  r, amp_cov,  t, sig=FALSE, pr = FALSE, ...) {
  
  
  n <- length(r)
  m <- sapply(r, length)
  
  sq <- logdet <- 0
  for (i in 1:n) {
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], param)
      
      U <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(U)) return(1e10) ## Exception handling;
    } else {
      
      S <- U <- diag(1, m[i])
    }
    rr <- as.numeric(r[[i]])
    
    sq <- sq +sum(backsolve(U, rr, transpose = TRUE)^2)
    logdet_tmp <- 0
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }
  
  sigmahat <- as.numeric(sq /sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  if (sig) return (sqrt(sigmahat))
  if (pr)  print(res)
  return(min(res, 1e10))
  
}


## Likehood using true Laplace approximation (as opposed to linearization)
likelihood.lap <- function(param, r, amp_cov, t, param.w, Zis, warp_cov, tw, sig=FALSE, pr = FALSE, w) {
  
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, param.w)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }
  
  n <- length(r)
  m <- sapply(r, length)
  
  sq <- logdet <- 0
  for (i in 1:n) {
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], param)
      
      U <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(U)) return(1e10) ## Exception handling;
    } 
    else {
      S <- U <- diag(1, m[i])
    }
    rr <- as.numeric(r[[i]])
    ZZ <- Zis[[i]]
    
    if (!is.null(warp_cov)) {
      A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
      LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
    } else {
      LR <- 0
    }
    sq <- sq + sum(backsolve(U, rr, transpose = TRUE)^2) 
    
    logdet_tmp <- 0
    if (!is.null(warp_cov)) logdet_tmp <- determinant(LR)$modulus[1]
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }
  
  if (!is.null(warp_cov)) {
    logdet <- logdet - n * determinant(Cinv)$modulus[1]
    sq <- sq + sum(backsolve(chol(C), w, transpose = TRUE)^2)
  } 
  
  sigmahat <- as.numeric(sq /sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  if (pr)  print(res)
  if (sig) return (sqrt(sigmahat))
  return(min(res, 1e10))
  
}

