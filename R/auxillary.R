




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


#' Spline weights 
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
#' @param design design matrix (as a list)
#' @param parallel Calculate the sum elements in parallel?
#' @rdname Spline_weights
#' 
#' @details 
#' Increasing splines might possibly not work together with designs.
#'
#' @return A matrix with spline weights
#' @export
#'
spline_weights <- function(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL, K, design, parallel = FALSE) {
  n <- length(y)
  if(length(design) != n) stop('Error! Unequal lengths of design and y')
  m <- sapply(y, nrow)
  des <- length(design[[1]])
  nb <- attr(basis_fct, 'df')*K*des
  
  if (is.null(weights)) weights <- rep(1, n)


    if (parallel) {
    dmat.dvec <-
      foreach(i = 1:n, yi = y, Sinvi = Sinv, ti = t, .combine = '+') %dopar% {

        basis <- diag(K) %x%  basis_fct(warp_fct(w[, i], ti))
        bSinv <- weights[i] * (t(basis) %*% Sinvi)
        c((design[[i]] %o% design[[i]]) %x% (bSinv %*% basis), design[[i]] %x% (bSinv %*% as.numeric(yi)))
      }
    Dmat <- matrix(dmat.dvec[1:(nb*nb)] , nb, nb)
    dvec <- matrix(dmat.dvec[length(Dmat)+1:nb], nb , 1)
  }
  else {
    Dmat <- matrix(0, nb, nb)
    dvec <- matrix(0, nb, 1)
    for (i in 1:n) {
      # Rewrite LS to reduce memory and calculation time.
      basis <- diag(K) %x%  basis_fct(warp_fct(w[, i], t[[i]]))
      bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
      Dmat <- Dmat + (design[[i]] %o% design[[i]]) %x% (bSinv %*% basis)
      dvec <- dvec + design[[i]] %x% (bSinv %*% as.numeric(y[[i]]))
    }
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
  else  c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  
  ce <- matrix(c, nc = K*des)
  return(ce)
}

splw.d <- spline_weights

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


## 'A' functions for exponential families

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



# c er matrix. Ens med men hurtigere end posterior2

#' Posterior likelihood 
#' 
#' @description Calculates the posterior likelihood for a single curve.
#'
#' @param w Warp values
#' @param warp_fct Warp function
#' @param t time points
#' @param y observations
#' @param basis_fct 
#' @param c matrix of spline coefficients
#' @param Sinv Precision matrix for amplitude
#' @param Cinv Precision matrix for w
#'
#' @return Value of posterior likelihood
#' @export
#'
posterior.lik <- function(w, warp_fct, t, y, basis_fct, c, Sinv, Cinv) {
  vt <- warp_fct(w, t)
  basis <- basis_fct(vt)
  r <- as.numeric(y - basis %*% c)
  return((t(r) %*% Sinv %*% r + t(w) %*% Cinv %*% w)[1])
}


