


#' Test
#'
#'@export
#'
ppMulti.test <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                       use.nlm = c(FALSE, FALSE), functional = NULL, 
                       amp_cov_par=NULL, paramMax = rep(T,length(amp_cov_par)),  warp_opt = TRUE, parallel.lik = FALSE,
                       like_optim_control = list(), pr=TRUE, design = NULL, inner_parallel = c(TRUE, TRUE), spammed = T,
                       save_temp = NULL, w0 = NULL, c0 = NULL) {
  
  require("spam")
  
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  em.iter <- if(is.na(iter[3])) 1 else iter[3]
  
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, nrow)
  K <- ncol(y[[1]])
  if (ncol(y[[1]]) == 1) warning("Use pavpop for one-dimensional functional objects")
  homeomorphisms <- 'soft'
  
  if (is.null(save_temp)) gem.tmp <- F
  else {
    gem.tmp <- T
    if (!is.character(save_temp)) stop("save_temp must be either False or a specified file location")
  }
  if(is.null(design)) stop("Not implemented without designs")
  
  ## S??rligt for denne
  if(parallel.lik) {
    print("Using parallelized likelihood")
    likelihood <- like.par
    if(!is.null(attr(amp_cov, "chol"))) print("Bruger smart choleski")
  }
  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'smooth') homeomorphisms <- 'no'
  
  # Unknown parameters
  warp_cov_par <- eval(attr(warp_cov, 'param'))
  n_par_warp <- length(warp_cov_par)
  n_par_amp <- length(amp_cov_par)
  
  
  p_warp <- if (!is.null(warp_cov) && warp_opt) 1:n_par_warp else c()
  
  ## Check if no. of ( lower) parameter limits correspond to ...
  
  if (!is.null(like_optim_control$lower) && length(like_optim_control$lower) > 1 && length(like_optim_control$lower) != n_par_amp + n_par_warp)
    warning("Mismatch between number of parameters and number of limits supplied! Problems may occur")
  
  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (!all(sapply(t, length) == m)) stop("Observations in y and t must have same length.")
  
  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]][,1])
    y[[i]] <- y[[i]][!missing_indices,]
    t[[i]] <- t[[i]][!missing_indices]
  }
  # Stored warped time
  t_warped <- t
  
  # Update m with cleaned data
  m <- sapply(y, nrow)
  
  
  
  if(!is.null(like_optim_control$parallel)) attr(amp_cov, "parallelLik") <- "T"
  # Initialize warp parameters
  if (is.null(w0)) w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  else w <- w0
  
  # Build amplitude covariances and inverse covariances
  if (is.null(amp_cov)) amp_cov <- diag_covariance
  
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))
  
  S <- Sinv <- S.chol <- list()
  
  S <- foreach (i = 1:n, ti = t, .noexport = c("t", "y")) %dopar% {
    Si <- amp_cov(ti, amp_cov_par)
    attr(Si, 'chol') <- chol(Si)
    attr(Si, 'inv') <- chol2inv(attr(Si, 'chol'))
    Si
  }
  
  for (i in 1:n) {
    S.chol[[i]] <- attr(S[[i]], 'chol')
    Sinv[[i]] <- attr(S[[i]], 'inv')
  }
  
  
  ## At gøre: Kombinér S.chol med inv_amp om muligt
  
  # Build warp covariance and inverse
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- solve(C)
  } else {
    C <- Cinv <- matrix(0, mw, mw)
  }
  
  
  # First estimate of  spline weights
  
  cis <- list() ## For designs  
  if (!is.null(design)) {
    c <- (splw.d(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
  } else {
    c <- splw(y, t, warp_fct, w, Sinv, basis_fct, K = K)
  }
  if (is.null(design)) cis[[i]] <- c
  else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
  
  # Construct warp derivative
  dwarp <- list()
  if (warp_type != 'smooth') {
    for (i in 1:n) {
      dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
      if (warp_type == 'piecewise linear') dwarp[[i]] <- as(dwarp[[i]], "dgCMatrix")
    }
  }
  
  # r and Z lists
  Zis <- list()
  r <- y
  
  # Initialize best parameters
  like_best <- Inf
  w_best <- w
  c_best <- c
  amp_cov_par_best <- amp_cov_par
  warp_cov_par_best <- warp_cov_par
  
  cat('Outer\t:\tInner \t:\tEstimates\n')
  for (iouter in 1:nouter) {
    if (halt_iteration & iouter != nouter) next
    # Outer loop
    if (iouter != nouter) cat(iouter, '\t:\t')
    
    
    wis <- list()
    
    for (iinner in 1:ninner) { ## Kriterium kan tilf??jes
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')
      
      
      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      w_res <- list()
      
      
      if (inner_parallel[1])  w_res <- ## Parallell prediction of warping parameters
        foreach(i = 1:n, Sinvi = Sinv, yi = y, tid = t, .noexport = c("Sinv", "S", "S.chol", "y", "t", "r", "Zis", "cis", "dwarp" )) %dopar% {

          ci <- if (!is.null(design)) c %*%  (design[[i]] %x% diag(K)) else c
          
          warp_optim_method <- 'CG'
          if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = tid, y = yi, c = ci, Sinv = Sinvi, Cinv = Cinv, basis_fct = basis_fct)$estimate 
          else  ww <- optim(par = w[, i], fn = posterior.lik, gr = NULL, method = warp_optim_method, warp_fct = warp_fct, t = tid, y = yi, c = ci, Sinv = Sinvi, Cinv = Cinv, basis_fct = basis_fct)$par
          
          if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
          return(ww)
        }
      else for ( i in 1:n) { ## Prediction as a for loop
        
        #if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
        #else cis[[i]] <- c 
        ci <- if (!is.null(design)) c %*%  (design[[i]] %x% diag(K)) else c
        
        warp_optim_method <- 'Nelder-Mead'
        if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = ci, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate 
        else  ww <- optim(par = w[, i], fn = posterior.lik, gr = NULL, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = ci, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
        
        if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
        w_res[[i]] <- ww
        
        if (spammed) eval(RZ.ting.spammed) 
        else eval(RZ.ting) 
      }
      
      ## Calculate R and Z   
      
      if (inner_parallel[1]) for (i in 1:n) {if (spammed) eval(RZ.ting.spammed) 
        else eval(RZ.ting) }
      
      ## Em algorithm stuff
      for (i.em in 1:em.iter) {
        
        if (i.em != em.iter) 
          for (i in 1:n)  {if (spammed) eval(RZ.ting.spammed) 
        else eval(RZ.ting) }
        
        # Update spline weights
        
        no.c <- length(c)
        sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
        west <- w
        
        # EM step
        cat(" EM algorithm! ")
        
        
        EM.expr <- expression({
          
          Amat <- matrix(0, no.c, no.c)
          
          rr <- as.numeric(rr)
          
          Sk <- backsolve(Scholi, ZZ, transpose = TRUE)
          Sz <-  t(Sk) %*% Sk
          
          west[,i] <- solve(Cinv +  Sz, t(Sk) %*% backsolve(Scholi, rr, transpose = TRUE))
          wvar <- sigma^2 * (C - C %*% (Sz - Sz %*% solve(Cinv +  Sz ,   Sz)) %*% C)
          

          R <- as.spam(t(desi) %x% diag(K)) %x% as.spam(bf(warpt))
          
          bd <-  bf( warpt, T)
          Rischol <- list()
          
          for (k in 1:mw) {
            Rk <- as.spam(t(desi) %x% diag(K)) %x% as.spam(wd[,k] *  bd )
            R <- R + Rk * (west[k,i] - w[k,i])
            
            Rischol[[k]] <- backsolve(Scholi, Rk, transpose = TRUE)
          }
          
          if (mw > 1) {
            for (l1 in 1:(mw-1)) {
              tr00 <- as.spam( (0.5*wvar[l1, l1]) *  Rischol[[l1]])
              for (l2 in (l1+1):mw) tr00 <- tr00 + wvar[l1, l2] *  Rischol[[l2]]
              Rs <- t(Rischol[[l1]]) %*% tr00
              Amat <- Amat + (Rs + t(Rs))
            }
            Amat <- Amat + wvar[mw,mw]*  (t(Rischol[[mw]]) %*%  Rischol[[mw]])
          } else stop("This method not implemented for mw < 2")
          
          SR <- backsolve(Scholi, R, transpose = TRUE)
          Amat <- Amat + t(SR) %*% SR
          bvek <- t(SR) %*% backsolve(Scholi, as.numeric(yi), transpose = TRUE)
          cbind(Amat, bvek)
        })
        
        ## Do the EM stuff
        
        if (inner_parallel[2]) Amatt <- 
          foreach (i = 1:n, rr = r, ZZ = Zis, Scholi = S.chol, wd = dwarp, yi= y, warpt = t_warped, desi = design, .combine = '+', 
                   .noexport = c("S", "Sinv", "Zis", "r", "Schol", "y", "t", "dwarp", "t_warped", "design"), 
                   .inorder = FALSE) %dopar% eval(EM.expr)
        else { ## for-loop
          Amatt <- matrix(0, no.c, no.c + 1)
          for (i in 1:n) {
            rr <- as.numeric(r[[i]])
            ZZ <- Zis[[i]]
            Scholi <- S.chol[[i]]
            wd <- dwarp[[i]]
            warpt <- t_warped[[i]]
            desi <- design[[i]]
            yi <- y[[i]]
            Amatt <- Amatt + eval(EM.expr)
          }
        }
        c.ny <- solve(as.spam(Amatt[, 1:no.c]), Amatt[, no.c + 1])
        
        if (pr) print(matrix(c.ny, nr = nrow(c)))
        c <- matrix(c.ny, nr = nrow(c), nc = ncol(c))
      }
      
      for (i in 1:n) {
        warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
        warp_change[2] <- max(warp_change[1], abs(w[, i] -  w_res[[i]]))
        w[, i] <- w_res[[i]]
      }
      if (warp_change[2] < 1e-2 / sqrt(mw)) break #TODO: Consider other criteria
      
    }
    for (i in 1:n)  {if (spammed) eval(RZ.ting.spammed) 
      else eval(RZ.ting) }
    
    # Likelihood estimation of parameters (outer loop)
    
    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      t_like <- t
      # if (warped_amp) t_like <- t_warped
      # Likelihood function
      par1 <- which(paramMax)
      parw <- n_par_amp + 1:n_par_warp
      
      like_eps <- if (is.null(like_optim_control$eps)) 1e-5 else like_optim_control$eps
      randomCycle <- if (is.null(like_optim_control$randomCycle)) -1 else like_optim_control$randomCycle
      optim.rule <- if (is.null(like_optim_control$optim.rule)) NULL else like_optim_control$optim.rule
      if (!is.null(optim.rule) && is.null(attr(optim.rule, "min"))) attr(optim.rule, "min") <- 0
      if (randomCycle[1] > -1 && randomCycle[2] < iouter) {
        par1 <- sample(par1, randomCycle[1])
        #        print(paste("Using random cycles. Optimizing on parameters",par1))
        print("Using random cycles. Optimizing on parameters ")
        cat(par1, "\n")
      } else if (!is.null(optim.rule) && attr(optim.rule, "min") < iouter) {
        op <- iouter %% length(optim.rule)
        if (op == 0) op <- length(optim.rule)
        par1 <- optim.rule[[op]]
        cat("Using cyclic optimization. Optimizing on parameters ",par1, "\n" )
      }
      
      if (n_par_warp > 0) like_fct <- function(pars) {
        
        par <- amp_cov_par
        if (warp_opt) {
          param.w <- pars[p_warp]
          par[par1]<- pars[- (p_warp)]
        }
        else {
          param.w <- warp_cov_par
          par[par1]<- pars
        }
        likelihood(par, param.w, r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw, pr = pr)
        
      } else like_fct <- function(pars) {
        par <- amp_cov_par
        par[par1]<- pars
        like.S(par, r = r, amp_cov = amp_cov, t = t, pr = pr)
      }
      
      # Likelihood gradient
      like_gr <- function(par) {
        
        res <- foreach(ip = 1:length(par), .combine = 'c') %:%
          foreach(sign = c(1, -1), .combine= '-') %dopar% {
            h <- rep(0, length(par))
            h[ip] <- sign * like_eps
            return(like_fct(par + h) / (2 * like_eps))
          }
        return(res)
      }
      
      # Estimate parameters using locally linearized likelihood
      lower  <- if (is.null(like_optim_control$lower)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$lower
      upper  <- if (is.null(like_optim_control$upper)) rep(Inf, n_par_amp + n_par_warp) else like_optim_control$upper
      method <- if (is.null(like_optim_control$method)) "L-BFGS-B" else like_optim_control$method
      ndeps <- if (is.null(like_optim_control$ndeps)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$ndeps
      maxit <- if (is.null(like_optim_control$maxit)) 20 else like_optim_control$maxit
      
      upper0  <- upper[c(1:n_par_warp , n_par_warp + par1)]
      lower0  <- lower[c(1:n_par_warp , n_par_warp + par1)]
      
      # Optim??r!
      cat("Optimization (outer loop) \n")
      paras <-  c(warp_cov_par, amp_cov_par[par1])
      #print(paras)
      
      if (use.nlm[1]) {  ## nlm optimization
        steptol <- if (is.null(like_optim_control$steptol)) 1e-6 else like_optim_control$steptol
        like_optim <- nlm.bound.xx(fct = like_fct , p = paras, lower = lower0, upper = upper0, init = TRUE, symmetric = TRUE, iterlim = maxit)
        param <- like_optim$estimate
        like_optim$value <- like_optim$minimum
      }
      else {    ## optim optimization
        like_optim <- optim(par = paras, like_fct, gr = like_gr, method = method, lower = lower0, upper = upper0, control = list(ndeps = ndeps, maxit = maxit))
        param <- like_optim$par
      }
      
      
      cat("Parameter values: ")
      
      if (!is.null(warp_cov) && warp_opt) warp_cov_par <- param[p_warp]
      if (!is.null(amp_cov)) amp_cov_par[par1] <- if (length(p_warp) > 0) param[- (p_warp)] else param
      
      if (like_optim$value <= like_best) {
        # Save parameters
        like_best <- like_optim$value
        w_best <- w
        c_best <- c
        amp_cov_par_best <- amp_cov_par
        warp_cov_par_best <- warp_cov_par
        
        cat(':\t', param, '\n')
        cat('Linearized likelihood:\t', like_best, '\n')
        if (gem.tmp) {
          cat('Saving estimates to ',save_temp, '\n')
          tmp_res = list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, sigma = 
                           likelihood(amp_cov_par_best, warp_cov_par_best, r, Zis, amp_cov, warp_cov, t, tw, sig=T),
                         warp_cov_par = warp_cov_par_best, like= like_best, iteration = iouter)
          save(tmp_res, file = save_temp)
        }
        
        # Update covariances
        S <- foreach (i = 1:n, ti = t, .noexport = c("t", "y")) %dopar% {
          Si <- amp_cov(ti, amp_cov_par)
          attr(Si, 'chol') <- chol(Si)
          attr(Si, 'inv') <- chol2inv(attr(Si, 'chol'))
          Si
        }
        for (i in 1:n) {
          S.chol[[i]] <- attr(S[[i]], 'chol')
          Sinv[[i]] <- attr(S[[i]], 'inv')
        }
        
        if (!is.null(warp_cov)) {
          C <- warp_cov(tw, warp_cov_par)
          Cinv <- solve(C)
        } else {
          C <- Cinv <- matrix(0, mw, mw)
        }
        
      } else {
        cat(':\tLikelihood not improved, returning best likelihood estimates.\n')
        halt_iteration <- TRUE
      }
    } else {
      # TODO: Should in principle be done before warps are updated in the final iteration!
      # Estimate of sigma if final iteration is reached
      if (nouter == 1) {
        w_best <- w
        c_best <- c
      }
      sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
      likeval <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=F)
      if (likeval <= like_best) like_best <- likeval
      else warning("Likelihood increased in final loop!")
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best))
}