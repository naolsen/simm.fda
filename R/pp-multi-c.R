

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
#' This is like ppMulti, but uses steps (default one step) of the EM algorithm to update spline coefficients in the linearized model.
#' More correct reults than ppMulti, but slower and differences are ususally small.
#' Requires spam package
#'
#' @param y List of observations in matrix form. NAs allowed
#' @param t List of corresponding time points. NAs not allowed
#' @param basis_fct Basis function for spline
#' @param warp_fct Warp function
#' @param amp_cov Amplitude covariance. Must be on the form function(t, param)
#' @param warp_cov Warp covariance 
#' @param iter two or three dimensional integer of maximal number of outer iterations &
#' maximal number of inner iterations per outer iteration and optionally maximal number of em updates (default = 1).
#' @param amp_cov_par Starting values for amplitude covariance parameters. There are no defaults.
#' @param warp_opt If FALSE, warp covariance parameters are kept fixed. 
#' @param paramMax Logical vector. Which amplitude parameters to optimise over? Defaults to all parameters.
#' @param parallel.lik Calculate likelihoods in parallel?
#' @param like_optim_control List of control options for optimization in outer loop. See details
#' @param use.nlm Use \code{nlm} instead of \code{optim} for optimization? First index for outer loop, second index for inner loop.
#' @param pr Printing option.
#' @param design May not be null
#' @param inner_parallel Should optimization of warps and matrices for EM algorithm be done in parallel?
#' @param save_temp Save estimates after each outer iteration? NULL or the file path.
#' @param w0 Starting values for warp. Should only be used if you have results from a previous run.
#' 
#' @details There has been less check on this function, so I cannot guarantee that it will behave as well as ppMulti.
#' Requires \code{spam} package, used for faster calculations for sparse matrices.
#' Regarding parallellization of EM inner step (matrix stuff): 
#' Results are combined in order they are returned by worker threads. Due to numericalities, this might imply that two runs are not sure to return completely identical results.
#' Also note that parallelization doesn't always work on Windows. (Use Linux if possible)
#'
#' @return A list of estimates
#' @export
#' @seealso \link{ppMulti}
#'
#' @examples See \link{ppMulti}
ppMulti.em <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                    use.nlm = c(FALSE, FALSE), functional = NULL, 
                    amp_cov_par=NULL, paramMax = rep(T,length(amp_cov_par)),  warp_opt = TRUE, parallel.lik = FALSE,
                    like_optim_control = list(), pr=TRUE, design = NULL, inner_parallel = c(TRUE, TRUE),
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
    #yvek[[i]] <- unlist(y[[i]])
    #tvek[[i]] <- rep(t[[i]], K)
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
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    
    S[[i]] <- amp_cov(t[[i]], amp_cov_par)
    S.chol[[i]] <- chol(S[[i]])
    
    if (inv_amp) {
      Sinv[[i]] <- inv_amp_cov(t[[i]], amp_cov_par)
    } else {
      Sinv[[i]] <- chol2inv(S.chol[[i]])
    }
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
          foreach(i = 1:n, Sinvi = Sinv, yi = y, tid = t, .noexport = c("Sinv", "S", "S.chol", "y", "t", "r", "Zis", "cis", "dwarp")) %dopar% {
            #if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
            #else cis[[i]] <- c 
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
          
        eval(RZ.ting) 
      }
      
      ## Calculate R and Z   
      
      if (inner_parallel[1]) for (i in 1:n)  eval(RZ.ting) 
      
      ## Em algorithm stuff
      for (i.em in 1:em.iter) {
      
      if (i.em != em.iter) 
        for (i in 1:n) eval(RZ.ting) 
        
      # Update spline weights

      no.c <- length(c)
      sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
      west <- w
      
      # EM step
      cat(" EM algorithm! ")
      
      
      EM.expr <- expression({
        
        Amat <- matrix(0, no.c, no.c)
        mi <- m[i]*K
        
        #Scholi <- S.chol[[i]]
        rr <- as.numeric(rr)
        #ZZ <- Zis[[i]]
        
        Sk <- backsolve(Scholi, ZZ, transpose = TRUE)
        Sz <-  t(Sk) %*% Sk
        
        west[,i] <- solve(Cinv +  Sz, t(Sk) %*% backsolve(Scholi, rr, transpose = TRUE))
        wvar <- sigma^2 * (C - C %*% (Sz - Sz %*% solve(Cinv +  Sz ,   Sz)) %*% C)
        
        #warpt <- t_warped[[i]]
        #Ri <- list()
        R <- as.spam(t(design[[i]]) %x% diag(K)) %x% as.spam(bf(warpt))
        
        #wd <- dwarp[[i]]
        bd <-  bf( warpt, T)
        Rischol <- list()
        
        for (k in 1:mw) {
          Rk <- as.spam(t(design[[i]]) %x% diag(K)) %x% as.spam(wd[,k] *  bd )
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
        foreach (i = 1:n, rr = r, ZZ = Zis, Scholi = S.chol, wd = dwarp, yi= y, warpt = t_warped, .combine = '+', 
                 .noexport = c("S", "Sinv", "Zis", "r", "Schol", "y", "t", "dwarp", "t_warped"), .inorder = FALSE) %dopar% eval(EM.expr)
      else { ## for-loop
        Amatt <- matrix(0, no.c, no.c + 1)
        for (i in 1:n) {
          rr <- as.numeric(r[[i]])
          ZZ <- Zis[[i]]
          Scholi <- S.chol[[i]]
          wd <- dwarp[[i]]
          warpt <- t_warped[[i]]
          yi <- y[[i]]
          Amatt <- Amatt + eval(EM.expr)
        }
      }
      c.ny <- solve(as.spam(Amatt[, 1:no.c]), Amatt[, no.c + 1])
      
     #  c.ny <- solve(Amat, bvek) Ikke ekspressivt
     # print(c.ny)
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
    for (i in 1:n)  eval(RZ.ting) 
    
    # Likelihood estimation of parameters (outer loop)
    
    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
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
          epsilon <- like_eps
          
          res <- rep(0, length(par))
          for (ip in  1:length(par)) {
            for (sign in c(1, -1)) {
              h <- rep(0, length(par))
              h[ip] <- sign * like_eps
              res[ip] <- res[ip] + sign * like_fct(par + h) / (2 * like_eps)
            }
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
          for (i in 1:n) {
            twarped <- t[[i]]
            
            S[[i]] <- amp_cov(twarped, amp_cov_par)
            S.chol[[i]] <- chol(S[[i]])
            
            if (inv_amp) {
              Sinv[[i]] <- inv_amp_cov(twarped, amp_cov_par)
            } else {
              Sinv[[i]] <- chol2inv(S.chol[[i]])
            }
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






#' Returns linearization stuff
#'
#' @param y List of observations in matrix form. NAs allowed
#' @param t List of corresponding time points. NAs not allowed
#' @param basis_fct Basis function for spline
#' @param warp_fct Warp function
#' @param amp_cov Amplitude covariance. Must be on the form function(t, param)
#' @param warp_cov Warp covariance 
#' @param iter Number of iterations on warp optimization. 0 allowed.
#' @param amp_cov_par Values for amplitude covariance parameters. 
#' @param w0 Starting values for warp. NULL is allowed.
#' @param pr Printing option
#' @param design Design for the experiments. Should be given as a list of one-dimensional vectors or as a design matrix.
#' @param inner_parallel Should the inner optimization be done parallelly?
#' @param eval_likelihood Calculate linearized likelihood of final estimate? Not yet implemented.
#'
#' @export
#' @return A list of results including r and z
#'
#' @seealso \link{ppMulti}
pavpop_returner_zogr <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = 5,
                       use.nlm = FALSE, functional = NULL, amp_cov_par=NULL,  w0 = NULL,  
                       pr=T, design = NULL, inner_parallel = FALSE, eval_likelihood = FALSE) {
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[1]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, nrow)
  K <- ncol(y[[1]])
  
  homeomorphisms <- 'soft'

  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'piecewise linear' & warp_type != 'smooth') homeomorphisms <- 'no'
  if (warp_type == 'identity') warp_cov <- NULL
  
  # Unknown parameters
  warp_cov_par <- eval(attr(warp_cov, 'param'))
  
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
  # Update m with cleaned data
  m <- sapply(y, nrow)
  
  cis <- list() ## For designs
  # Design part. If matrix, convert to list
  if (is.matrix(design)) {
    des <- design
    design <- list()
    for (i in 1:n) design[[i]] <- des[i,]
  }
  if (!is.null(design) && length(design) != n) stop("design must have same length or number of rows as the length of y.")
  
  ## Check if functionality package is used
  
  if (!is.null(functional)) {
    require("fctbases")
    
    bf0 <- functional()
    basis_fct <- setup.functional.basis(bf0, t[[1]][1], isTRUE(formals(functional)$ownDeriv))
    
    on.exit({ ## Close connection after use.
      print("Closing")
      fctbases::removeMember(bf0)
    })
  }
  
  # Stored warped time
  t_warped <- t
  
  
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
  like_best <- NULL
  sigma <- NULL
  w_best <- w
  c_best <- c
  amp_cov_par_best <- amp_cov_par
  warp_cov_par_best <- warp_cov_par
  
  cat('Inner \t:\tEstimates\n')
    
    
    wis <- list()
    
    if (ninner > 0)     #  wis <- list()
      for (iinner in 1:ninner) {
        # Inner loop
        cat(iinner, '\t')
        
        # Predict warping parameters for all functional samples
        warp_change <- c(0, 0)
        if (homeomorphisms == 'hard') {
          #TODO: constrainOptim
          stop("Hard homeomorphic constrained optimization for warps is not implemented.")
        } 
        else if (attr(warp_fct, "mw") != 0) {
          # Parallel prediction of warping parameters
          gr <- NULL
          w_res <- list()
          
          if (inner_parallel)  w_res <- 
            foreach(i = 1:n) %dopar% {
              if (!is.null(functional))  stop("functional && inner_parallel does not work together!")

              warp_optim_method <- 'CG'
              if (use.nlm) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = y[[i]], 
                                     c = if (!is.null(design)) c %*%  (design[[i]] %x% diag(K)) else c, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate 
              else  ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], 
                                c = if (!is.null(design)) c %*%  (design[[i]] %x% diag(K)) else c, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
              
              if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
              return(ww)
            }
          else for ( i in 1:n) {
            
            if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
            else cis[[i]] <- c 
            warp_optim_method <- 'Nelder-Mead'
            if (use.nlm) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate 
            else  ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
            
            if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
            w_res[[i]] <- ww
          }
          
          for (i in 1:n) {
            warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
            warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
            w[, i] <- w_res[[i]]
          }
        }
        else cat(" Skipping inner optimization")
        # Update spline weights
        if (is.null(design)) {
          c <- splw(y, t, warp_fct, w, Sinv, basis_fct, K = K)
        } else {
          c <- (splw.d(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
        }
        
        
        
        if (warp_change[2] < 1e-2 / sqrt(mw)) break 
        
      }
    
    ### Outer loop part. 
    
    ## Construct residual vector for given warp prediction
    
    Zis <- list()
    r <- y
    
    for (i in 1:n)  eval(RZ.ting) 
    
    w_best <- w
    c_best <- c
    
    if(eval_likelihood) {
      sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
      like_best <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=F)
    } 
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best, Zis = Zis, r = r))
}


