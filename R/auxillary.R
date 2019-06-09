




RZ.ting <- expression({
  if (is.null(design)) cis[[i]] <- c
  else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
  
  # Compute warped time
  twarped <- t_warped[[i]] <- warp_fct(w[, i], t[[i]])
  if (!is.null(warp_cov)) {
    if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
    
    Zis[[i]] <- multi.Zi(twarped, dwarp[[i]], basis_fct, cis[[i]], mw) ## Opdateret. multi.Zi er hurtigere.
  } else {
    Zis[[i]] <- matrix(0, m[i]*K, mw)
  }
  
  ## Opdateret. Hurtigere evaluering af y - r - Zw^0
  if (nrow(w) != 1) {
    r[[i]] <- y[[i]] - as.numeric(basis_fct(twarped) %*% cis[[i]]) + as.numeric(Zis[[i]] %*% w[, i])
  }
  else {
    rrr <- y[[i]]
    
    for (k in 1:K) {
      rrr[,k] <- rrr[,k] - basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
    }
    r[[i]] <- rrr
  }
  
})

RZ.ting.spammed <- expression({
  if (is.null(design)) cis[[i]] <- c
  else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
  
  # Compute warped time
  twarped <- t_warped[[i]] <- warp_fct(w[, i], t[[i]])
  if (!is.null(warp_cov)) {
    if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
    
    Zis[[i]] <- as.spam(multi.Zi(twarped, dwarp[[i]], basis_fct, cis[[i]], mw)) ## Opdateret. multi.Zi er hurtigere.
  } else {
    Zis[[i]] <- matrix(0, m[i]*K, mw)
  }
  
  ## Opdateret. Hurtigere evaluering af y - r - Zw^0
  if (nrow(w) != 1) {
    r[[i]] <- y[[i]] - as.numeric(basis_fct(twarped) %*% cis[[i]]) + as.numeric(Zis[[i]] %*% w[, i])
  }
  else {
    rrr <- y[[i]]
    
    for (k in 1:K) {
      rrr[,k] <- rrr[,k] - basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
    }
    r[[i]] <- rrr
  }
  
})

#' Linearised likelihood function
#'
#' @param param Amplitude parameters
#' @param param.w Warp parameters
#' @param r,Zis Residual vectors/matrices
#' @param amp_cov Amplitude covariance function
#' @param warp_cov Warp covariance function
#' @param t Observation time points
#' @param tw 
#' @param sig Return sigma or likelihood?
#' @param pr Print option. Only use for debugging.
#' @rdname likelis
#' 
#' @details like.par performs parallelized likelihood
#'
#' @return -2*logL(parameters)
#' @export
#'
likelihood <- function(param, param.w, r, Zis, amp_cov, warp_cov, t, tw, sig=FALSE, pr = FALSE) {
  
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
      # U <- chol(S)
      # HANDLE ERRORS:
      U <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-5, m[i])))
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
like.par <- function (param, param.w, r, Zis, amp_cov, warp_cov, t, tw, sig=FALSE, pr = FALSE) {
  
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, param.w)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }
  
  n <- length(r)
  m <- sapply(r, length)
  
  sqLogd <- foreach(i = 1:n, ZZ = Zis, rr = r, tid = t, .combine = '+', .noexport = c("Zis", "r", "t"),
     .inorder = FALSE) %dopar% {
        if (!is.null(amp_cov)) {
          S <- amp_cov(tid, param)
          # U <- chol(S)
          # HANDLE ERRORS:
          U <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-5, m[i])))
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




## Spline weights for individual designs. 

#' Spline weights
#'
#' @rdname Spline_weights
#' 
#' @export
#'
splw.d <- function(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL, K, design = list()) {
  
  n <- length(y)
  if(length(design) != n) stop('Error! Unequal lengths of design and y')
  m <- sapply(y, nrow)
  des <- length(design[[1]])
  nb <- attr(basis_fct, 'df')*K*des
  
  if (is.null(weights)) weights <- rep(1, n)
  
  Dmat <- matrix(0, nb, nb)
  dvec <- matrix(0, nb, 1)
  
  for (i in 1:n) {
    basis <- t(design[[i]]) %x% diag(K) %x%  basis_fct(warp_fct(w[, i], t[[i]])) # I think this is the first time I have used that the kronecker product is associative
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% as.numeric(y[[i]])
  }
  
  c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  
  ce <- matrix(c, nc = K*des)
  return(ce)
}


## Splinev??gte til flere dimensioner med kryds-korrelation

#' Spline weights for individual designs. 
#' @description  Spline weights
#'
#' @param y list of observations
#' @param t list of time points
#' @param warp_fct warping function
#' @param w warp parameters
#' @param Sinv precision matricees
#' @param basis_fct Basis function
#' @param weights weights (optional)
#' @param K dimension
#' @param design design
#' @rdname Spline_weights
#' 
#' @details splw is used when no design is supplied.
#' splw.d is used when a desigen is supplied.
#' splw has been updated to handle increasing splines.
#'
#' @return A matrix with spline weights
#' @export
#'
splw <- function(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL, K) {
  n <- length(y)
  m <- sapply(y, nrow)
  nb <- attr(basis_fct, 'df')*K
  
  if (is.null(weights)) weights <- rep(1, n)
  
  Dmat <- matrix(0, nb, nb)
  dvec <- matrix(0, nb, 1)
  
  for (i in 1:n) {
    basis <-diag(K) %x%  basis_fct(warp_fct(w[, i], t[[i]]))
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% as.numeric(y[[i]])
  }
  
  if (attr(basis_fct, 'increasing')) {
    
    intercept <- !(attr(basis_fct, 'intercept')) 
    
    #  tryCatch({ ## To be improved...
    #    c <- solve.QP(Dmat, dvec, Amat = diag(x = c(intercept,  rep(1, nb-1 )) ))$solution
    #  }, error = function(e) {
    
    #    c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
    #    cej0 <- (c != 0)
    #    c[cej0] <- solve.QP(Dmat[cej0, cej0], dvec[cej0], Amat = diag(x = c(intercept,  rep(1, length(cej0) - 1) )))$solution
    
    #diag(Dmat) <- diag(Dmat) +  1e-6*mean(diag(Dmat))
    #c <- solve.QP(Dmat, dvec, Amat = diag(x = c(intercept,  rep(1, nb-1 )) ))$solution
    rank <- rankMatrix(Dmat)[1]
    if (rank != nb) {
      index <- nb
      indices <- 1:nb
      
      # Perhaps this can be done smarter or better?
      while (length(indices) > rank) {
        tmp_indices <- indices[indices != index]
        if (rankMatrix(Dmat[tmp_indices, tmp_indices]) == rank) {
          indices <- tmp_indices
        }
        index <- index - 1
      }
      c <- rep(0, ncol(basis))
      c[indices] <- solve.QP(Dmat = Dmat[indices, indices], dvec = dvec[indices,],
                             Amat = diag(nrow = length(indices)))$solution
      
    }
    else  c <- solve.QP(Dmat, dvec, Amat = diag(x = c(intercept,  rep(1, nb-1 )) ))$solution
    
    #  })
  }
  else c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  
  ce <- matrix(c, nc = K)
  return(ce)
}

## noget med nul
splw.new <- function(y, tw, arp_fct, w, Sinv, basis_fct, weights = NULL, K) {
  n <- length(y)
  m <- sapply(y, nrow)
  nb <- attr(basis_fct, 'df')*K
  
  if (is.null(weights)) weights <- rep(1, n)
  
  Dmat <- matrix(0, nb, nb)
  dvec <- matrix(0, nb, 1)
  
  for (i in 1:n) {
    basis <-diag(K) %x%  basis_fct(warp_fct(w[, i], t[[i]]))
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% as.numeric(y[[i]])
  }
  
  if (attr(basis_fct, 'increasing')) {
    
    intercept <- !(attr(basis_fct, 'intercept')) 
    diag(Dmat)[diag(Dmat) == 0] <- 1  ## 
    
    rank <- rankMatrix(Dmat)[1]
    if (rank != nb) {
      index <- nb
      indices <- 1:nb
      
      # Perhaps this can be done smarter or better?
      while (length(indices) > rank) {
        tmp_indices <- indices[indices != index]
        if (rankMatrix(Dmat[tmp_indices, tmp_indices]) == rank) {
          indices <- tmp_indices
        }
        index <- index - 1
      }
      c <- rep(0, ncol(basis))
      c[indices] <- solve.QP(Dmat = Dmat[indices, indices], dvec = dvec[indices,],
                             Amat = diag(nrow = length(indices)))$solution
      
    }
    else  c <- solve.QP(Dmat, dvec, Amat = diag(x = c(intercept,  rep(1, nb-1 )) ))$solution
    
  }
  else c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  
  ce <- matrix(c, nc = K)
  return(ce)
}


## 'A' functions

### Neg.binom

#' Create negative binomial
#'
#' @param r 
#'
#' @export
#'
#' @examples Afkt <- create.negbin(4.5)
#' #simfd.ed(..., ed.fct = Afkt, ...) ## Use with simfd.ed
#' 
#' @rdname Afkt
create.negbin <- function(r) {
  
  Afkt <- function(x,y) (r+y)*log1p(exp(x)/r) ## Consistent with Poisson model; i.e. r \pil infty gives Poisson.
  attr(Afkt, 'diff2') <-  function(x, y) (y+r) * exp(x)/(r*(1+exp(x)/r)^2)
  
  return(Afkt)
}  

#' Create Poisson model
#'
#' @export
#'
#' @examples Afkt <- create.poisson()
#' #simfd.ed(..., ed.fct = Afkt, ...) ## Use with simfd.ed
#' @rdname Afkt
create.poisson <- function() {
  
  Afkt <- function(x,y) exp(x)
  attr(Afkt, 'diff2') <-  function(x, y) exp(x)
  
  return(Afkt)
} 


