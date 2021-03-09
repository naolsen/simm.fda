





#' Likelihood in ED models
#'
#' @param param Amplitude parameters
#' @param param.w Warp covariance parameters
#' @param r,Zis Residual vectors/matrices
#' @param A.diff2 2nd derivative for exponential family
#' @param amp_cov 
#' @param warp_cov 
#' @param t Observation time points
#' @param tw 
#' @param pr Print option
#'
#' @return A number
#' 
likelihood.ed <- function (param, param.w, r, Zis, A.diff2, amp_cov, warp_cov, t, tw, pr = FALSE) {

  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, param.w)
    Cinv <- chol2inv(chol(C))
  }
  else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }
  n <- length(r)
  m <- sapply(r, length)
  sq <- logdet <- 0
  for (i in 1:n) {
    
    rr <- as.numeric(r[[i]])
    ZZ <- Zis[[i]]

    S <- amp_cov(t[[i]], param)
    if (!is.null(warp_cov)) {
      S <- S + ZZ %*% C %*% Matrix::t(ZZ)
    }
    U <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(U)) return(1e10) ## Exception handling;

    sq <- sq + sum(backsolve(U, rr, transpose = TRUE)^2)
    
    
    ## logdet(S) + logdet(tilde{\Sigma}) = logdet(I + SD)
    S <- 2*A.diff2[[i]] * S + diag(m[i])
    
    logdet <- logdet + determinant(S)$modulus
    
  }
  as.numeric(sq + logdet)
}


posterior.lik_u <- function(u, y, vt, Sinv, A.fct) { ## alpha = 1/r
  0.5*t(u- vt) %*% Sinv %*% (u-vt) +
    sum(A.fct(u, y) - y*u)
}


#' Simultaneous inference for misaligned functional data in an exponential family setting
#' 
#' This function applies the simm.fda methodology in a setting of certain exponential models (e.g. poisson) followiong the approach outlined in article[..]
#'
#' @param y List of observations. NAs allowed
#' @param t List of corresponding time points. NAs not allowed
#' @param basis_fct Basis function for spline
#' @param warp_fct Warp function
#' @param ed_fct Function defining the exponential family. Must have attribute diff2. See \link{Afkt} for examples.
#' @param amp_cov Amplitude covariance function. Must be on the form \code{function(t, param)}
#' @param warp_cov Warp covariance function. Must be on the form \code{function(t, param)}
#' @param iter two-dimensional integer of maximal number of outer iterations &
#' maximal number of inner iterations per outer iteration.
#' @param w0 Starting values for predicted warps. Should only be used if you have results from a previous run.
#' @param u0 Starting values for predicted trajectories. Should only be used if you have results from a previous run. 
#' @param use.nlm se \code{nlm} instead of \code{optim} for optimization? First index for outer loop, second index for inner loop.
#' @param suppressLik Suppress if likelihood has increased
#' @param amp_cov_par Starting values for amplitude covariance parameters. There are no defaults.
#' @param paramMax Logical vector. Which amplitude parameters to optimise over? Defaults to all parameters.
#' May be overwritten by supplying control parameters.
#' @param warp_opt If FALSE, warp covariance parameters are kept fixed. 
#' @param like_optim_control  List of control options for optimization in outer loop. See \link{ppMulti} for details.
#' @param pr Printing option.
#' @param design Design for the experiments. Should be given as a list of one-dimensional vectors or as a design matrix.
#' @param inner_parallel  Should the inner optimization be done parallelly?
#' @param save_temp  Save estimates after each outer iteration? NULL or the file path.
#' 
#' @details sim.fd can only handle one-dimensional functional data. It should NOT be applied to Gaussian data, unless sigma^2 is known and fixed beforehand.
#' Not all relevant control parameters have been checked for compliance; however parallelization does work (on Linux)
#'
#' 
#'
#' @return A list of estimates and predictions of w and u
#' @export
#'
#' @seealso \link{ppMulti}
#'
simfd.ed <- ppMulti.ed <- function(y, t, basis_fct, warp_fct, ed_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                          w0 = NULL, u0 = NULL, use.nlm = c(FALSE, FALSE), suppressLik = FALSE,
                          amp_cov_par=NULL, paramMax = rep(TRUE, length(amp_cov_par)), warp_opt = TRUE,
                          like_optim_control = list(), design = NULL, inner_parallel = FALSE, save_temp = NULL) {
  
  ## Ops�tnings-ting
  
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, length)
  
  if (is.null(save_temp)) gem.tmp <- F
  else {
    gem.tmp <- T
    if (!is.character(save_temp)) stop("save_temp must be either NULL or a specified file location")
  }
  
  Ad2 <- attr(ed_fct, 'diff2')
  
  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type == 'identity') warp_cov <- NULL
  
  
  # Unknown parameters
  if (!is.null(warp_cov)) {
    warp_cov_par <- eval(attr(warp_cov, 'param'))
  }
  else {
    warp_cov_par <- c()
    mw <- length(tw)
  }
  n_par_warp <- length(warp_cov_par)
  
  if (mw == 0)  {
    n_par_warp <- 0
    warp_opt <- FALSE
  }
  
  p_warp <- if (!is.null(warp_cov) && warp_opt) 1:n_par_warp else c()
  
  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (!all(sapply(t, length) == m)) stop("Observations in y and t must have same length.")
  
  
  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]])
    y[[i]] <- y[[i]][!missing_indices]
    t[[i]] <- t[[i]][!missing_indices]
  }
  
  t_warped <- t  # Stored warped time
  m <- sapply(y, length) # Update m with cleaned data
  
  # Design part. If matrix, convert to list
  if (is.null(design))
    design <- as.list(rep(1, n))
  else if (is.matrix(design)) {
    des <- design
    design <- list()
    for (i in 1:nrow(des)) design[[i]] <- des[i,] ## Slightly nicer with nrow
  }
  if (length(design) != n) stop("design must have same length or number of rows as the length of y.")
  cis <- list()
  
  # Build amplitude covariances and inverse covariances
  if (is.null(amp_cov)) {
    amp_cov <- function(t, par) diag(length(t))
    amp_cov_par <- c()
    paramMax <- logical(0)
  }
  n_par_amp <- length(amp_cov_par)
  
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))
  
  S <- Sinv <- list()
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    
    S[[i]] <- amp_cov(t[[i]], amp_cov_par)
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
  
  # Initialize warp parameters
  if (is.null(w0))   w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  else w <- w0
  if (is.null(like_optim_control$warp_crit)) like_optim_control$warp_crit <- 1e-2
  if (is.null(like_optim_control$optim_w)) like_optim_control$optim_w <- TRUE
  
  vts <- list() 
  
  Adiffs <- list()
  Avals <- rep(0, n)
  
  # Initialize  posterior values
  u <- list()
  for (i in 1:n) { ## Kan parallelliseres
    
    u[[i]] <- if (!is.null(u0)) u0[[i]] else
      optim(rep(0, m[i]), posterior.lik_u, y = y[[i]], Sinv = Sinv[[i]], vt = rep(0, m[i]), A.fct = ed_fct)$par
    Avals[i] <- sum(ed_fct(u[[i]], y[[i]])- y[[i]]*u[[i]])
    Adiffs[[i]] <- Ad2(u[[i]], y[[i]])
    
  }

  # Check if no. of (lower) parameter limits correspond to no. of parameters
  if (!is.null(like_optim_control$lower) && length(like_optim_control$lower) > 1 && length(like_optim_control$lower) != n_par_amp + n_par_warp)
    warning("Mismatch between number of parameters and number of limits supplied! Problems may occur")
  
  # Estimate spline weights
  c <- spline_weights(u, t, warp_fct, w, Sinv, basis_fct, K = 1, design=design)
  for (i in 1:n) cis[[i]] <- c %*% design[[i]]
  
  # Construct warp derivative
  dwarp <- list()
  if (warp_type != 'smooth') {
    for (i in 1:n) {
      dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
      if (warp_type == 'piecewise linear') dwarp[[i]] <- as(dwarp[[i]], "dgCMatrix")
    }
  }
  r <- Zis <- list()
  
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
    
    
    #  wis <- list()
    for (iinner in 1:ninner) {
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')
      
      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      
      if (attr(warp_fct, "mw") != 0) {
        # Parallel prediction of warping parameters
        gr <- NULL
        w_res <- list()
        if (like_optim_control$optim_w) {
          
          if (inner_parallel)  w_res <- 
            foreach(i = 1:n, Sinvi = Sinv, yi = y, ui = u, .noexport = c("y", "Sinv", "S", "dwarp", "r", "Zis", "cis", "dwarp")) %dopar% {

              ci <- if (!is.null(design)) c %*% design[[i]]  else c

              if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = ui, c = ci, Sinv = Sinvi, Cinv = Cinv, basis_fct = basis_fct)$estimate
              else  ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = 'CG', warp_fct = warp_fct, t = t[[i]], y = ui, c = ci, Sinv = Sinvi, Cinv = Cinv, basis_fct = basis_fct)$par

              return(ww)
            }
          else for ( i in 1:n) {
            cis[[i]] <- c %*% design[[i]]

            if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = u[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate
            else  ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = 'Nelder-Mead', warp_fct = warp_fct, t = t[[i]], y = u[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
            w_res[[i]] <- ww
          }
          
          for (i in 1:n) {
            warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
            warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
            w[, i] <- w_res[[i]]
          }}
        
        for (i in 1:n) {
          vts[[i]] <- basis_fct(warp_fct(w[,i], t[[i]])) %*% cis[[i]]
          
          u[[i]] <- optim(u[[i]], posterior.lik_u, y = y[[i]], Sinv = Sinv[[i]], vt = vts[[i]], A.fct = ed_fct)$par
          Avals[i] <- sum(ed_fct(u[[i]], y[[i]])- y[[i]]*u[[i]])
          Adiffs[[i]] <- Ad2(u[[i]], y[[i]])
        }
        
        
      }
      else cat("No warping function provided! Skipping inner optimization.")
      # Update spline weights
      c <- spline_weights(u, t, warp_fct, w, Sinv, basis_fct, K = 1, design=design)
      
      if (warp_change[2] < like_optim_control$warp_crit / sqrt(mw) && like_optim_control$optim_w) break 
      
    }
    
    ### Outer loop part. 
    
    ## 1. Construct residual vector for given warp prediction
    
    Zis <- r <- list()

    for (i in 1:n) {
      cis[[i]] <- c %*% (design[[i]])
      twarped <- t_warped[[i]] <- warp_fct(w[, i], t[[i]])
      if (!is.null(warp_cov)) {
        if (warp_type == "smooth") 
          dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        Zis[[i]] <- simm.fda:::multi.Zi(twarped, dwarp[[i]], basis_fct, 
                                        cis[[i]], mw)
      }
      else {
        Zis[[i]] <- matrix(0, m[i], mw)
      }

      if (nrow(w) != 1) {
        r[[i]] <- u[[i]] - as.numeric(basis_fct(twarped) %*% 
                                        cis[[i]]) + as.numeric(Zis[[i]] %*% w[, i])
      }
      else {
        rrr <- u[[i]]
        for (k in 1:1) {
          rrr[, k] <- rrr[, k] - basis_fct(twarped) %*% cis[[i]][, k] + Zis[[i]][(m[i] * (k - 1) + 1):(m[i] * k), ] * w[, i]
        }
        r[[i]] <- rrr
      }
    }
    
    ##for (i in 1:n) eval(substitute(RZ.ting, list(y = u)) ?? G?r noget!!
    
    
    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      
      # Likelihood function
      par1 <- which(paramMax)
      
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
      
      
      Ad2vals <- list()
      for (i in 1:n) Ad2vals[[i]] <- Ad2(u[[i]], y[[i]])
      
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
        likelihood.ed(par, param.w, r = r, Zis = Zis, A.diff2 = Ad2vals, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw)
        
      }
      else stop("not implemented without warp optimization!")
      
      
      # Likelihood gradient
      like_gr <- function(par) {
        
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
      
      ## Estimate parameters using locally linearized likelihood
      # Control parameters:
      lower  <- if (is.null(like_optim_control$lower)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$lower
      upper  <- if (is.null(like_optim_control$upper)) rep(Inf, n_par_amp + n_par_warp) else like_optim_control$upper
      method <- if (is.null(like_optim_control$method)) "L-BFGS-B" else like_optim_control$method
      ndeps <- if (is.null(like_optim_control$ndeps)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$ndeps
      maxit <- if (is.null(like_optim_control$maxit)) 20 else like_optim_control$maxit
      
      upper0  <- upper[c(p_warp , n_par_warp + par1)]
      lower0  <- lower[c(p_warp , n_par_warp + par1)]
      
      # Optim??r!
      cat("Optimization (outer loop) \n")
      paras <- if (warp_opt) c(warp_cov_par, amp_cov_par[par1]) else amp_cov_par[par1]
      
      
      if (use.nlm[1]) {  ## nlm optimization
        steptol <- if (is.null(like_optim_control$steptol)) 1e-6 else like_optim_control$steptol
        like_optim <- nlm.bound(fct = like_fct , p = paras, lower = lower0, upper = upper0, iterlim = maxit)
        param <- like_optim$estimate
        like_optim$value <- like_optim$minimum
      }
      else {    ## optim optimization
        like_optim <- optim(par = paras, like_fct, gr = like_gr, method = method, lower = lower0, upper = upper0, control = list(ndeps = ndeps, maxit = maxit))
        param <- like_optim$par
      }
      
      like_optim$value <- like_optim$value + 2*sum(Avals) ## Tilf?j log-lik fra ed-familie.
      
      cat("Parameter values: ")
      
      if (!is.null(warp_cov) && warp_opt) warp_cov_par <- param[p_warp]
      if (!is.null(amp_cov)) amp_cov_par[par1] <- if (length(p_warp) > 0) param[- (p_warp)] else param
      
      if (like_optim$value <= like_best || suppressLik) {
        # Save parameters
        like_best <- like_optim$value
        w_best <- w
        c_best <- c
        amp_cov_par_best <- amp_cov_par
        warp_cov_par_best <- warp_cov_par
        itermax <- iouter
        
        cat('\t', param, '\n')
        cat('Linearized likelihood:\t', like_best, '\n')
        if (gem.tmp) {
          cat('Saving estimates to ',save_temp, '\n')
          tmp_res <- list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best,
                         warp_cov_par = warp_cov_par_best, like= like_best, iteration = iouter, u = u)
          save(tmp_res, file = save_temp)
        }
        # Update covariances
        for (i in 1:n) {
          twarped <- t[[i]]
          
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
        
      } else {
        itermax <- iouter - 1
        cat(':\tLikelihood not improved, returning best likelihood estimates.\n')
        if (!suppressLik) halt_iteration <- TRUE
      }
    } else {
      # Estimate of sigma if final iteration is reached
      if (nouter == 1) {
        w_best <- w
        c_best <- c
      }
      #sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, u = u, 
              like = like_best, iterations = itermax))
}



## NOTE: different parametrization
## Denne skal v?re i brug

##Tilføjes senere.
#likelihood.ed2 
#ppMulti.ed2     





