

## 

## Parametre:
## paramMax: Hvilke af amplitudeparametrene skal der maksimeres over?  upper & lower skal passe med disse.
## w.c.p: pt ikke i brug. Skal v?re TRUE.

## win: s?t til falsk for ikke at parallelisere warp pr?diktioner. Kan v?re n?dvendigt for at undg? crash
## Parallel registering ikke laengere en del af ppMulti. Skal goeres inden kald af ppMulti.

## save_temp: Gem undervejs? Pas p? med parallellisering
## Vigtig ??ndring: Nu skal gr??nser angives korrekt ift. alle parametre, dvs length(upper) = n_par_warp + n_par_amp. De ikke-brugte gr??nser m?? gerne v??re tomme
## Tilsvarende for lower (faktisk nok mere intuitivt). aendringen er kun relevant hvis man ikke maksimerer over alle parametre.


#' Simulataneous inference for Misaligned Multivarariate functional data
#'
#' @param y list of observations in matrix form. NAs allowed
#' @param t list of corresponding time points
#' @param basis_fct Basis fct.
#' @param warp_fct Warp fct.
#' @param amp_cov Amplitude covariance. Must be on the form function(t, param)
#' @param warp_cov Warp covariance 
#' @param iter two-dimensional integer of maximal number ofouter iterations 
#' and maxiamal number of inner iterations per outer iteration.
#' @param parallel 
#' @param amp_cov_par Starting values for amplitude covariance parameters. There are no defaults.
#' @param paramMax logical vecter. Which amplitude parameters to optimise over? Defaults to all parameters.
#' @param w.c.p Do not change
#' @param like.alt Do not change either
#' @param like_optim_control List of control options for optimization in outer loop. See details
#' @param pr print option
#' @param design Design for the experiments. Should be given as a list of one-dimensional vectors.
#' @param win 
#' @param save_temp Save estimates after each outer iteration? FALSE, NULL or the file path.
#'
#' @details ppMulti returns error if applied on one-dimensional functional data.
#' 
#' Control parameters
#' lower, upper, method, ndeps, maxit will be sent to optim. See optim for more details. 
#' The first entries of lower/upper correspond to warp parameters, while the rest corresponds to
#' amplitude parameters. ppMulti does match upper/lower with corresponding entries in amp_cov, which is important when not all parameters are maximized over. 
#' This is for consistency with randomCycle and optimRule.
#' 
#' randomCycle and optimRule are two ways of optimizing on only a subset of the parameters at a time. TBD: descriptions of these.
#'
#' @aliases simm.fda
#'
#' @seealso 
#' Important details can be found in simm-fda-short-desc.pdf
#'
#' @return A list of estimates
#' @export
#'
#' @examples \link{example}
#' 
#' \donttest{
#' # Data originates from http://mocap.cs.cmu.edu/
#'
#' ## Data is provided as
#' # MCD.data y-values
#' # MCD.time time points
#
# # The time points have been scaled to be within [0.1, 0.9] ...
#'
#' # Make basis function
#' bf <- make_basis_fct(seq(0, 1, len=32)[2:31], intercept=T, boundary = c(0, 1))
#'
#' ## Make warp function
#' tw <- seq(0, 1, length = 5)[2:4] ## anchor points for hyman spline
#' wf.noshift <- make_warp_fct(type="smooth", tw=tw)
#'
#' ## warp function with shift
#' wf <- w.shift(wf.noshift, 0.25)
#'
#' warp_cov <- warp_and_shift_cov(c(0.5, 1))
#'
#'
#' ## For this example we assume no amplitude cross-correlation, otherwise ...
#' # Matern covariance; smoothness parameter 2, unknown range.
#'
#' wrapper <- function(t, par) {
#'   mvMatern(t, par[1], 2, rep(1,3), diag(par[2:4]))
#' }
#'
#' lower <- c(1e-4, 1e-4,rep(1e-4, 4))
#' upper <- rep(1e5, 6)
#'
#' # Note this will take a long time:
#' mcd.res <-ppMulti(MCD.data, MCD.time, bf, wf, amp_cov = wrapper, warp_cov, amp_cov_par = c(0.2, rep(100, 3)), pr = F, paramMax = rep(T, 4),
#'                   like_optim_control = list(lower = lower, upper = upper), win = F, iter = c(10, 10))
#'}
#'

ppMulti <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
                    parallel = list(n_cores = 1, parallel_likelihood = FALSE), 
                    amp_cov_par=NULL, paramMax = rep(T,length(amp_cov_par)),  w.c.p = T, like.alt = F,
                    like_optim_control = list(), pr=F, design = NULL, win =F, save_temp = NULL) {
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, nrow)
  K <- ncol(y[[1]])
  if (ncol(y[[1]]) == 1) warning("simm.fda cannot be expected to work on one-dimensional curves.")
  homeomorphisms <- 'soft'
  if (is.null(save_temp)) gem.tmp <- F
  else {
    gem.tmp <- T
    if (!is.character(save_temp)) stop("save_temp must be either False or a specified file location")
  }
  
  cis <- list() ## For designs
  ## problemer med parallelitet
  if (win) {
    tryFejl <- FALSE
  } else tryFejl <- TRUE
  
  if(like.alt) {
    print("Skifter funktion")
    likelihood <- like.par
    if(!is.null(attr(amp_cov, "chol"))) print("Bruger smart choleski")
  }
  
  
  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'piecewise linear' & warp_type != 'smooth') homeomorphisms <- 'no'
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
  
  
  if (mw == 0)  n_par_warp <- 0
  n_par_amp <- length(amp_cov_par)
  
  p_warp <- if (!is.null(warp_cov)) 1:n_par_warp else c()
  
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

  
  # Build amplitude covariances and inverse covariances
  if (is.null(amp_cov)) amp_cov <- diag_covariance
  
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
  w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  
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
    for (iinner in 1:ninner) {
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')
      
      
      
      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      if (homeomorphisms == 'hard') {
        #TODO: constrainOptim
        stop("Hard homeomorphic constrained optimization for warps is not implemented.")
      } 
      else if (attr(warp_fct, "mw") != 0) {
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
              return(ww)
            }
            tryPar <- F
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
            

            if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
            w_res[[i]] <- ww
            wis[[i]] <- w_res[[i]]
          }
          
          tryFejl <- T
        }
        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
        }
        
        for (i in 1:n) w[, i] <- w_res[[i]]
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
    
    
    # Likelihood estimation of parameters (outer loop)
    # Pre-compute warped time
    for (i in 1:n) {
      t_warped[[i]] <- warp_fct(w[, i], t[[i]])
    }
    
    # If the amplitude variation is assumed to be varying in warped time, the amplitude covariances are updated
    
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
        Zis[[i]] <- matrix(0, m[i]*K, mw)
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
    
    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      t_like <- t

      
      # Likelihood function
      par1 <- which(paramMax)
      parw <- n_par_amp + p_warp
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
                par[par1]<- pars[- (p_warp)]
        if (w.c.p) param.w <- pars[p_warp]
        else param.w <- warp_cov_par
        likelihood(par, param.w, r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw, pr = pr)
        
      }
      else like_fct <- function(pars) {
        par <- amp_cov_par
        par[par1]<- pars
        like.S(par, r = r, amp_cov = amp_cov, t = t, pr = pr)
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
      
      upper0  <- upper[c(p_warp , n_par_warp + par1)]
      lower0  <- lower[c(p_warp , n_par_warp + par1)]
      
      # Optim??r!
      cat("Optimization (outer loop) \n")
      paras <-  c(warp_cov_par, amp_cov_par[par1])
      
      like_optim <- optim(par = paras, like_fct, gr = like_gr, method = method, lower = lower0, upper = upper0, control = list(ndeps = ndeps, maxit = 20))
      param <- like_optim$par
  #    print(param)
      cat("Parameter values: ")
      
      if (!is.null(warp_cov)) warp_cov_par <- param[p_warp]
      if (!is.null(amp_cov)) 
      amp_cov_par[par1] <- if (n_par_warp > 0) param[- (p_warp)] else param
      
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
        cat('\t', param, '\n')
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
## randomCycle: vektor af laengde to. [1] angiver antal tilfaeldigt udtrukne paramertre;
## [2] hvornaar cyklen skal begynde
## optim.rule: cyclisk optimering. liste bestaaende af hvilke parametre der skal opdateres i hvilken blok
## Alle parametre boer vaere indeholdt og der maa gerne vaere gentagelser. attr(, 'min') kan bruges til at saette hvornaar opim.rule skal begynde
## Bemaerk at de to ovenstaaende overskriver paramMax .

# c er matrix. Ens med men hurtigere end posterior2

#' Posterior likelihood 
#' 
#' @description Calculates the posterior likelihood for a single curve.
#'
#' @param w warp values
#' @param warp_fct 
#' @param t time points
#' @param y observations
#' @param basis_fct 
#' @param c matrix of spline coefficients
#' @param Sinv Precision matrix for amplitude
#' @param Cinv Precision matrix for w
#'
#' @return
#' @export
#'
posterior.lik <- function(w, warp_fct, t, y, basis_fct, c, Sinv, Cinv) {
  vt <- warp_fct(w, t)
  basis <- basis_fct(vt)
  r <- as.numeric(y - basis %*% c)
  return((t(r) %*% Sinv %*% r + t(w) %*% Cinv %*% w)[1])
}

## Splinev??gte til flere dimensioner med kryds-korrelation

## Spline weights for individual designs. 
#' @description  Spline weights
#'
#' @param y list of observations
#' @param t list of time points
#' @param warp_fct 
#' @param w 
#' @param Sinv precision matricees
#' @param basis_fct 
#' @param weights weights (optional)
#' @param K dimension
#' @param design design
#' @rdname Spline_weights
#' 
#' @details splw is used when no design is supplied.
#' splw.d is used when a desigen is supplied.
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
    intercept <- attr(basis_fct, 'intercept')
    indices <- 1:nb
    rank <- rankMatrix(Dmat)[1]
    index <- 1
    
    # Perhaps this can be done smarter or better?
    while (length(indices) > rank) {
      tmp_indices <- indices[indices != index]
      if (rankMatrix(Dmat[tmp_indices, tmp_indices]) == rank) {
        indices <- tmp_indices
      }
      index <- index + 1
    }
    c <- rep(0, ncol(basis))
    c[indices] <- solve.QP(Dmat = Dmat[indices, indices],
                           dvec = dvec[indices,],
                           Amat = diag(nrow = length(indices)))$solution
  } else {
    c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  }
  ce <- matrix(c, nc = K)
  return(ce)
}

## Spline weights for individual designs. 
#' Spline weights
#'
#' @rdname Spline_weights
#' 
#' @return
#' @export
#'
#' @examples
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
    basis <- t(design[[i]]) %x% diag(K) %x%  basis_fct(warp_fct(w[, i], t[[i]])) # I think this is the first time I have used that the kronecker product is commutative
    
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% as.numeric(y[[i]])
  }
  
  c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  
  ce <- matrix(c, nc = K*des)
  return(ce)
}



#' Linearised likelihood function
#'
#' @param param 
#' @param param.w 
#' @param r 
#' @param Zis 
#' @param amp_cov 
#' @param warp_cov 
#' @param t 
#' @param tw 
#' @param sig 
#' @param pr Print value?
#'
#' @return -2*logL(parameters)
#' @export
#'
#' @examples
likelihood <- function(param, param.w, r, Zis, amp_cov, warp_cov, t, tw, sig=F, pr = F) {

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
  if (sig) return (sqrt(sigmahat))
  if (pr)  print(res)
  return(min(res, 1e10))
  
}
## Used when no warp is given
like.S <- function(param,  r, amp_cov,  t,  sig=F, pr = F) {
  

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
    
    sq <- sq +sum(backsolve(U, rr, transpose = TRUE)^2)
    logdet_tmp <- 0
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }
  if (!is.null(warp_cov)) logdet <- logdet 
  
  sigmahat <- as.numeric(sq /sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  if (sig) return (sqrt(sigmahat))
  if (pr)  print(res)
  return(min(res, 1e10))
  
}


## Allow for smart decompositions. "Non-specified" covariances not supported
like.par <- function (param, param.w, r, Zis, amp_cov, warp_cov, t, tw, sig=F, pr = T) {
  
  
  C <- warp_cov(tw, param.w)
  Cinv <- chol2inv(chol(C))
  
  n <- length(r)
  m <- sapply(r, length)
  
  chFunk <- attr(amp_cov, "chol")
  parall <- !is.null(attr(amp_cov, "parallelLik"))
  
  if(is.null(chFunk)) chSmart <- F
  else chSmart <- T
  
  
  n <- length(r)
  m <- sapply(r, length)
  
  if(parall) {
    registerDoParallel(cores = n)
  }
  
  
  sqLogd <-   foreach (i = 1:n) %dopar% {
    
    
    if (chSmart) U <- chFunk(t[[i]], param)
    else {
      S <- amp_cov(t[[i]], param)
      # HANDLE ERRORS:
      U <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-6, m[i])))
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
    sqq <-  (sum(backsolve(U, rr, transpose = TRUE)^2)
             - t(x) %*% LR %*% x)
    logdet_tmp <- 0
    if (!is.null(warp_cov)) logdet_tmp <- determinant(LR)$modulus[1]
    logd <-  - (logdet_tmp - 2 * sum(log(diag(U))))
    c(sqq, logd)
  }
  if (!is.null(warp_cov)) logdetC <- - n * determinant(Cinv)$modulus[1]
  # print(sqLogd)
  
  sqLogd <- matrix(unlist(sqLogd), nr=2)

  sigmahat <- as.numeric(sum(sqLogd[1,]) /sum(m))
  #res <- sum(m) * log(sigmahat) + logdet
  res <- sum(m) * log(sigmahat) + sum(sqLogd[2,]) + logdetC
  if (sig) return (sqrt(sigmahat))
  if (pr)  print(res)
  return(min(res, 1e10))
  
}   



