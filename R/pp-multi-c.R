

## 

## Parametre:
## paramMax: Hvilke af amplitudeparametrene skal der maksimeres over?  upper & lower skal passe med disse.
## w.c.p: pt ikke i brug. Skal v?re TRUE.

## pr: printoption
## design: Hvis der bruges design, e.g. en kurve for hver person. Skal v?re samme design p? alle koordinater.
## win: s?t til falsk for ikke at parallelisere warp pr?diktioner. Kan v?re n?dvendigt for at undg? crash

## save_temp: Gem undervejs? Pas p? med parallellisering
## Vigtig ??ndring: Nu skal gr??nser angives korrekt ift. alle parametre, dvs length(upper) = n_par_warp + n_par_amp. De ikke-brugte gr??nser m?? gerne v??re tomme
## Tilsvarende for lower (faktisk nok mere intuitivt). ?ndringen er kun relevant hvis man ikke maksimerer over alle parametre.

#' ppMulti with update steps using EM algorithm
#' 
#' This is like ppMulti, but uses one step of the EM algorithm to update spline coefficients in the linearized model.
#' Theoretically better than ppMulti, but slower and the differences should be neglible.
#' Requires spam package
#'
#' @param y 
#' @param t 
#' @param basis_fct 
#' @param warp_fct 
#' @param amp_cov 
#' @param warp_cov 
#' @param iter 
#' @param parallel 
#' @param amp_cov_par 
#' @param paramMax 
#' @param w.c.p 
#' @param like.alt 
#' @param like_optim_control 
#' @param pr 
#' @param design May not be null
#' @param win 
#' @param inner_control 
#' @param save_temp 
#' @param w0 
#' 
#' @details There has been less check on this function, so I cannot guarantee that it will behave as well as ppMulti.
#' Requires \code{spam} package, used for faster calculations for sparse matrices.
#'
#' @return
#' @export
#' @seealso \link{ppMulti}
#'
#' @examples
ppMulti.em <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                    parallel = list(n_cores = 1, parallel_likelihood = FALSE),
                    amp_cov_par=NULL, paramMax = rep(T,length(amp_cov_par)),  w.c.p = T, like.alt = F,
                    like_optim_control = list(), pr=T, design = NULL, win =T, inner_control = list(),
                    save_temp = NULL, w0 = NULL) {
  
  require("spam")
  
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, nrow)
  K <- ncol(y[[1]])
  if (ncol(y[[1]]) == 1) stop("Use pavpop for one-dimensional functional objects")
  homeomorphisms <- 'soft'
  if (is.null(save_temp)) gem.tmp <- F
  else {
    gem.tmp <- T
    if (!is.character(save_temp)) stop("save_temp must be either False or a specified file location")
  }
  if(is.null(design)) stop("Not implemented without designs")
  
  cis <- list() ## For designs
  ## problemer med parallelitet
  if (win) {
    tryFejl <- FALSE
  } else tryFejl <- TRUE
  ## S??rligt for denne
  
  
  if(like.alt) {
    print("Skifter funktion")
    likelihood <- like9
    if(!is.null(attr(amp_cov, "chol"))) print("Bruger smart choleski")
  }
  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'piecewise linear' & warp_type != 'smooth') homeomorphisms <- 'no'
  
  # Unknown parameters
  warp_cov_par <- eval(attr(warp_cov, 'param'))
  n_par_warp <- length(warp_cov_par)
  n_par_amp <- length(amp_cov_par)
  
  ## Check if no. of ( lower) parameter limits correspond to ...
  
  if (!is.null(like_optim_control$lower) && length(like_optim_control$lower) > 1 && length(like_optim_control$lower) != n_par_amp + n_par_warp)
    warning("Mismatch between number of parameters and number of limits supplied! Problems may occur")
  
  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (!all(sapply(t, length) == m)) stop("Observations in y and t must have same length.")
  
  yvek <- list()
  tvek <- list()
  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]][,1])
    y[[i]] <- y[[i]][!missing_indices,]
    t[[i]] <- t[[i]][!missing_indices]
    yvek[[i]] <- unlist(y[[i]])
    tvek[[i]] <- rep(t[[i]], K)
  }
  # Stored warped time
  t_warped <- t
  
  # Update m with cleaned data
  m <- sapply(y, nrow)
  
  # Initialize cluster
  if (win) registerDoParallel(cores = parallel$n_cores)
  
  if(!is.null(like_optim_control$parallel)) attr(amp_cov, "parallelLik") <- "T"
  # Initialize warp parameters
  if (is.null(w0)) w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  else w <- w0
  
  # Build amplitude covariances and inverse covariances
  if (is.null(amp_cov)) amp_cov <- diag_covariance
  
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))
  
  S <- Sinv <- list()
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    
    S[[i]] <- amp_cov(t[[i]], amp_cov_par)
    #print(dim(S[[i]]))
    if (inv_amp) {
      Sinv[[i]] <- inv_amp_cov(t[[i]], amp_cov_par)
    } else {
      Sinv[[i]] <- chol2inv(chol(S[[i]]))
    }
  }
  
  # Build warp covariance and inverse
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- solve(C)
  } else {
    C <- Cinv <- matrix(0, mw, mw)
  }
  
  # First estimate of  spline weights
  
  if (!is.null(design)) {
    c <- (splw.d(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
  } else {
    c <- splw(y, t, warp_fct, w, Sinv, basis_fct, K = K)
  }
  
  # Construct warp derivative
  dwarp <- list()
  if (warp_type != 'smooth') {
    for (i in 1:n) {
      dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
      if (warp_type == 'piecewise linear') dwarp[[i]] <- as(dwarp[[i]], "dgCMatrix")
    }
  }
  
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
    
    if (is.null(inner_control$warp_change_crit)) inner_control$warp_change_crit <- T

    
    #posteriorIncrease <- rep(T, n)      
    posteriorLik <- rep(Inf, n)
    #  pSum <- Inf
    
    wis <- list()
    
    for (iinner in 1:ninner) { if (T) { ## Kriterium kan tilf??jes
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')
      
      
      
      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      if (homeomorphisms == 'hard') {
        #TODO: constrainOptim
        stop("Hard homeomorphic constrained optimization for warps is not implemented.")
      } else {
        # Parallel prediction of warping parameters
        tryPar <- T  ## Pr?v parallelk?rsek
        w_res <- try (
          if (!tryFejl) {
            foreach(i = 1:n) %dopar% {
              gr <- NULL
              if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
              else cis[[i]] <- c 
              # warp_optim_method <- 'Nelder-Mead'
              warp_optim_method <- 'CG'
              ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
              if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
              tryPar <- F
              return(ww)
            }
          })
        if (tryPar) { ## .. ellers ..
          w_res <- wis
          gr <- NULL
          
          for ( i in 1:n) {
            if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
            else cis[[i]] <- c 
            warp_optim_method <- 'Nelder-Mead'
            ww0 <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)
            ww <- ww0$par
            posteriorLik <- ww0$val
            
            if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
            w_res[[i]] <- ww
            wis[[i]] <- w_res[[i]]
          }

          tryFejl <- T
        }
        
        
        
        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
          w[, i] <- w_res[[i]]
        }
        
      }
      
      Zis <- list()
      r <- y
      
      
      ## Her er det lidt mindre klart
      for (i in 1:n) {
        t_warped[[i]] <- warp_fct(w[, i], t[[i]])
        if (is.null(design)) cis[[i]] <- c
        else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
        twarped <- t_warped[[i]]
        if (!is.null(warp_cov)) {
          if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
          Zis[[i]] <- matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, 1]), m[i], mw)
          for (k in 2:K) {
            Zis[[i]] <- rbind(Zis[[i]], matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, k]), m[i], mw))
          }
          
        } else {
          Zis[[i]] <- Matrix(0, m[i]*K, mw)
        }
        
        
        rrr <- y[[i]]
        for (k in 1:K) {
          
          if (nrow(w) == 1) {
            rrr[,k] <- rrr[,k] - basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
          } else {
            rrr[,k] <- rrr[,k] -  basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] %*% w[, i]
          }
        }
        r[[i]] <- rrr
      }
      # Update spline weights

      no.c <- length(c)
      Amat <-  matrix(0, no.c, no.c)
      
      bvek <- rep(0, no.c)
      sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
      if (T) likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, pr = T)
      west <- w
      
      # EM algorithm step
      val0 <-c(0,0)
      for (i in 1:n) {
       if (pr) cat(" ", i)
        
        mi <- m[i]*K
         Sinvi <- Sinv[[i]]
        rr <- as.numeric(r[[i]])
        ZZ <- Zis[[i]]
        
         tZZSz <- t(ZZ) %*% Sinvi
         Sz <- tZZSz %*% ZZ
        west[,i] <- solve(Cinv +  Sz , tZZSz %*% rr )
       A <- ZZ %*% solve(Cinv +  Sz ,   tZZSz )

          wvar <- sigma^2 * (C - C %*% (Sz - tZZSz %*%  A %*% ZZ) %*% C)
          
          if (pr) print(wvar)
          if(pr) print (c(west[,i], w[,i]))
        
          warpt <- t_warped[[i]]
        Ri <- list()
        R <- as.spam(t(design[[i]]) %x% diag(K)) %x% as.spam(bf(warpt))
        
       # wd <-  warp_fct(w[,i], t[[i]], w_grad = T)
        wd <- dwarp[[i]]
        bd <-  bf( warpt, T)
        
        for (k in 1:mw) {
          Ri[[k]] <- as.spam(t(design[[i]]) %x% diag(K)) %x% as.spam(wd[,k] *  bd )
          R <- R +  Ri[[k]] * (west[k,i] - w[k,i])
        }
        
        if (F && i == 21) { R_test <- matrix(nr = mi, nc = mw)
        for (k in 1:mw) R_test[,k] <- Ri[[k]] %*% as.numeric(c)
        print(cbind(R_test, Zis[[i]]))
         for (k in 1:mw) R_test[,k] <-  as.spam(wd[,k] *  bd ) %*% c[,1:3]
          print(cbind(R_test, Zis[[i]]))
        #Zi(t, wd, basis_fct , )
        }
        
        for (l1 in 1:mw) {
          #Rs <- as.spam(t(Ri[,,l1])) %*% Sinvi
          Rs <- t(Ri[[l1]]) %*% Sinvi
          tr00 <- as.spam(matrix(0, mi , no.c))
          for (l2 in 1:mw) {
            tr00 <- tr00 + wvar[l1, l2] *  Ri[[l2]]# Ri[,,l2]
          
          }
        # tr_temp <- tr_temp + Rs %*% tr00
          Amat <- Amat +  Rs %*% tr00
        }

        
              ## R vil ikke anerkende det som symmetriske matricer
        bas <-  t(R) %*% Sinvi
        Amat <- Amat + bas %*% R 
        bvek <- bvek + bas %*% as.numeric(y[[i]])      
      }
      
      Amat <- as.spam(Amat)
      Amat <- 0.5 * (Amat + t(Amat) )   
         # print(mean(abs(Amatuden-Amat)))
      c.ny <- solve(Amat, bvek)
     # print(c.ny)
      if (pr) print(matrix(c.ny, nr = nrow(c)))
      c <- matrix(c.ny, nr = nrow(c), nc = ncol(c))
    }    
      if (warp_change[2] < 1e-2 / sqrt(mw)) break #TODO: Consider other criteria

    }
  
    
    
    # Likelihood estimation of parameters (outer loop)
    
    
    
    # Pre-compute warped time
    for (i in 1:n) {
      t_warped[[i]] <- warp_fct(w[, i], t[[i]])
    }
    
    
    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    
    
    ## Her er det lidt mindre klart
    for (i in 1:n) {
      if (is.null(design)) cis[[i]] <- c
      else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
      twarped <- t_warped[[i]]
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        Zis[[i]] <- matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, 1]), m[i], mw)
        for (k in 2:K) {
          Zis[[i]] <- rbind(Zis[[i]], matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, k]), m[i], mw))
        }
        
      } else {
        Zis[[i]] <- Matrix(0, m[i]*K, mw)
      }
      
      
      rrr <- y[[i]]
      for (k in 1:K) {
        
        if (nrow(w) == 1) {
          rrr[,k] <- rrr[,k] - basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
        } else {
          rrr[,k] <- rrr[,k] -  basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] %*% w[, i]
        }
      }
      r[[i]] <- rrr
    }
    
    ## If warps are strange, try to re-run inner loop and compare. warp_cov has to be.
    
    
    
    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      t_like <- t
      # if (warped_amp) t_like <- t_warped
      # Likelihood function
      par1 <- which(paramMax)
      parw <- n_par_amp + 1:n_par_warp
      #print(par1)
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
      
      like_fct <- function(pars) {
        
        par <- amp_cov_par
        
        par[par1]<- pars[- (1:n_par_warp)]
        if (w.c.p) param.w <- pars[1:n_par_warp]
        else param.w <- warp_cov_par
        likelihood(par, param.w, r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw, pr = pr)
        
      }
      
      # Likelihood gradient
      like_gr <- NULL
      if (parallel$parallel_likelihood) {
        # Construct parallel gradient
        like_gr <- function(par) {
          epsilon <- 1e-5
          rep(1:length(par), each = 2)
          res <- foreach(ip = 1:length(par), .combine = 'c') %:%
            foreach(sign = c(1, -1), .combine= '-') %dopar% {
              h <- rep(0, length(par))
              h[ip] <- sign * epsilon
              return(like_fct(par + h) / (2 * epsilon))
            }
          return(res)
        }
      } else {
        # Construct parallel gradient
        like_gr <- function(par) {
          epsilon <- 1e-5
          rep(1:length(par), each = 2)
          res <- rep(0, length(par))
          for (ip in  1:length(par)) {
            for (sign in c(1, -1)) {
              h <- rep(0, length(par))
              h[ip] <- sign * epsilon
              res[ip] <- res[ip] + sign * like_fct(par + h) / (2 * epsilon)
            }
          }
          return(res)
        }
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
      cat("Optimization (outer loop)" , parw, "\n")
      paras <-  c(warp_cov_par, amp_cov_par[par1])
      #print(paras)
      
      like_optim <- optim(par = paras, like_fct, gr = like_gr, method = method, lower = lower0, upper = upper0, control = list(ndeps = ndeps, maxit = 20))
      param <- like_optim$par
      print(param)
      
      if (!is.null(warp_cov)) warp_cov_par <- param[1:n_par_warp]
      if (!is.null(amp_cov)) amp_cov_par[par1] <- param[- (1:n_par_warp)]
      
      if (like_optim$value <= like_best) {
        # Save parameters
        like_best <- like_optim$value
        w_best <- w
        c_best <- c
        amp_cov_par_best <- amp_cov_par
        warp_cov_par_best <- warp_cov_par
        
        # Update covariances
        for (i in 1:n) {
          twarped <- t[[i]]
          #  if (warped_amp) twarped <- t_warped[[i]]
          # Check if an amplitude covariance is defined
          S[[i]] <- amp_cov(twarped, amp_cov_par)
          if (inv_amp) {
            Sinv[[i]] <- inv_amp_cov(twarped, amp_cov_par)
          } else {
            Sinv[[i]] <- chol2inv(chol(S[[i]]))
          }
        }
        
        if (!is.null(warp_cov)) {
          C <- warp_cov(tw, warp_cov_par)
          Cinv <- solve(C)
        } else {
          C <- Cinv <- matrix(0, mw, mw)
        }
        cat(':\t', param, '\n')
        cat('Linearized likelihood:\t', like_best, '\n')
        if (gem.tmp) {
          cat('Saving estimates to ',save_temp, '\n')
          tmp_res = list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, sigma = 
                           likelihood(amp_cov_par_best, warp_cov_par_best, r, Zis, amp_cov, warp_cov, t, tw, sig=T),
                         warp_cov_par = warp_cov_par_best, like= like_best, iteration = iouter)
          save(tmp_res, file = save_temp)
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
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best))
}

## Kontrolparametre: alle defaulter til NULL
## maxit: maximum af iterationer i \tt{optim} i ydre loop
## randomCycle: vektor af l?ngde to. [1] angiver antal tilf?ldigt udtrukne paramertre;
## [2] hvorn?r cyklen skal begynde
## optim.rule: cyclisk optimering. liste best?ende af hvilke parametre der skal opdateres i hvilken blok
## Alle parametre b?r v?re indeholdt og der m? gerne v?re gentagelser. attr(, 'min') kan bruges til at s?tte hvorn?r opim.rule skal begynde
## Bem?rk at de to ovenst?ende overskriver paramMax .






pavpop_returner_zogr <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                       parallel = list(n_cores = 1, parallel_likelihood = FALSE), use_warp_gradient = FALSE,
                       amp_cov_par=NULL,  w0 = NULL,  pr=T, design = NULL, win =T, inner_control = list()) {
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, nrow)
  K <- ncol(y[[1]])
  if (ncol(y[[1]]) == 1) stop("Use pavpop for one-dimensional functional objects")
  homeomorphisms <- 'soft'

  
  
  cis <- list() ## For designs
  ## problemer med parallelitet
  if (win) {
    tryFejl <- FALSE
  } else tryFejl <- TRUE
  ## S??rligt for denne
  
  
  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'piecewise linear' & warp_type != 'smooth') homeomorphisms <- 'no'
  
  # Unknown parameters
  warp_cov_par <- eval(attr(warp_cov, 'param'))
  n_par_warp <- length(warp_cov_par)
  n_par_amp <- length(amp_cov_par)
  
  ## Check if no. of ( lower) parameter limits correspond to ...
  
  #if (!is.null(like_optim_control$lower) && length(like_optim_control$lower) > 1 && sum(paramMax)+ n_par_warp !=  length(like_optim_control$lower))
  #  print("Warning: there is a mismatch in number of (optimization) parameters and number of limits supplied!")
  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (!all(sapply(t, length) == m)) stop("Observations in y and t must have same length.")
  
  yvek <- list()
  tvek <- list()
  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]][,1])
    y[[i]] <- y[[i]][!missing_indices,]
    t[[i]] <- t[[i]][!missing_indices]
    yvek[[i]] <- unlist(y[[i]])
    tvek[[i]] <- rep(t[[i]], K)
  }
  # Stored warped time
  t_warped <- t
  
  # Update m with cleaned data
  m <- sapply(y, nrow)
  
  
  # Initialize warp parameters
  if (is.null(w0)) w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  else w <- w0
  
  # Build amplitude covariances and inverse covariances
  if (is.null(amp_cov)) amp_cov <- diag_covariance
  
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))
  
  S <- Sinv <- list()
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    
    S[[i]] <- amp_cov(t[[i]], amp_cov_par)
    #print(dim(S[[i]]))
    if (inv_amp) {
      Sinv[[i]] <- inv_amp_cov(t[[i]], amp_cov_par)
    } else {
      Sinv[[i]] <- chol2inv(chol(S[[i]]))
    }
  }
  
  # Build warp covariance and inverse
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- solve(C)
  } else {
    C <- Cinv <- matrix(0, mw, mw)
  }
  
  # Estimate spline weights
  
  if (!is.null(design)) {
    c <- (splw4(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
  } else {
    c <- splw3(y, t, warp_fct, w, Sinv, basis_fct, K = K)
  }
  
  # Construct warp derivative
  dwarp <- list()
  if (warp_type != 'smooth') {
    for (i in 1:n) {
      dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
      if (warp_type == 'piecewise linear') dwarp[[i]] <- as(dwarp[[i]], "dgCMatrix")
    }
  }
  
  # Initialize best parameters
  like_best <- Inf
  w_best <- w
  c_best <- c
  amp_cov_par_best <- amp_cov_par
  warp_cov_par_best <- warp_cov_par
  
  cat('Outer\t:\tInner \t:\tEstimates\n')
  for (iouter in 1:1) {
    if (halt_iteration & iouter != nouter) next
    # Outer loop
    if (iouter != nouter) cat(iouter, '\t:\t')
    
    
    wis <- list()
    
    for (iinner in 1:ninner) {
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')
      
      
      
      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      if (homeomorphisms == 'hard') {
        #TODO: constrainOptim
        stop("Hard homeomorphic constrained optimization for warps is not implemented.")
      } else {
        # Parallel prediction of warping parameters
        tryPar <- T  ## Pr?v parallelk?rsek
        w_res <- try (
          if (!tryFejl) {
            foreach(i = 1:n) %dopar% {
              gr <- NULL
              if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
              else cis[[i]] <- c 
              # warp_optim_method <- 'Nelder-Mead'
              warp_optim_method <- 'CG'
              ww <- optim(par = w[, i], fn = posterior3, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct, K=K)$par
              if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
              tryPar <- F
              return(ww)
            }
          })
        if (tryPar) { ## .. ellers ..
          w_res <- wis
          gr <- NULL
          
          for ( i in 1:n) {
            if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
            else cis[[i]] <- c 
            warp_optim_method <- 'Nelder-Mead'
            ww0 <- optim(par = w[, i], fn = posterior3, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct, K=K)
            ww <- ww0$par
            

            if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
            w_res[[i]] <- ww
            wis[[i]] <- w_res[[i]]
          }
          # print(sum(posteriorLik))
          # print(posteriorLik)
          
          tryFejl <- T
        }
        
        
        
        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
        }
        
        for (i in 1:n) w[, i] <- w_res[[i]]
      }
      
      # Update spline weights
      if (is.null(design)) {
        c <- splw3(y, t, warp_fct, w, Sinv, basis_fct, K = K)
      } else {
        c <- (splw4(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
      }
      
      
      
      if (inner_control$warp_change_crit && warp_change[2] < 3e-3 / sqrt(mw)) break #TODO: Consider other criteria
    }
    
    
    
    # Likelihood estimation of parameters (outer loop)
    
    
    
    # Pre-compute warped time
    for (i in 1:n) {
      t_warped[[i]] <- warp_fct(w[, i], t[[i]])
    }
    
    
    
    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    
    
    ## Her er det lidt mindre klart
    for (i in 1:n) {
      if (is.null(design)) cis[[i]] <- c
      else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
      twarped <- t_warped[[i]]
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        Zis[[i]] <- matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, 1]), m[i], mw)
        for (k in 2:K) {
          Zis[[i]] <- rbind(Zis[[i]], matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, k]), m[i], mw))
        }
        
      } else {
        Zis[[i]] <- Matrix(0, m[i]*K, mw)
      }
      
      
      rrr <- y[[i]]
      for (k in 1:K) {
        
        if (nrow(w) == 1) {
          rrr[,k] <- rrr[,k] - basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
        } else {
          rrr[,k] <- rrr[,k] -  basis_fct(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] %*% w[, i]
        }
      }
      r[[i]] <- rrr
    }
    
    ## If warps are strange, try to re-run inner loop and compare. warp_cov has to be.
    
    
    
    # Check wheter the final outer loop has been reached
    
    # TODO: Should in principle be done before warps are updated in the final iteration!
    # Estimate of sigma if final iteration is reached
    
    w_best <- w
    c_best <- c
    
    sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
  }
return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best, Zis = Zis, r = r))
}

#' simm.fda
#' 
#' @param t 
#' @param basis_fct 
#' @param warp_fct 
#' @param amp_cov 
#' @param warp_cov 
#' @param iter 
#' @param parallel 
#' @param amp_cov_par 
#' @param paramMax 
#' @param w.c.p 
#' @rdname simm.fda
#' 
#' 
#' 

