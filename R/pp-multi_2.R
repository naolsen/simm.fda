
# Som den skal v?re; kan overf?res til bibliotek.
## Recent changes: Added feature warp_opt. Some parallel stuff de-implemented. Small cleanups
# version 1.5
# 4-4-17: More cleanup and and small imprevements. Added argument inner_parallel

## Parametre:
## paramMax: Hvilke af amplitudeparametrene skal der maksimeres over?  upper & lower skal passe med disse.

## win: s?t til falsk for ikke at parallelisere warp pr?diktioner. Kan v?re n?dvendigt for at undg? crash
## Parallel registering ikke laengere en del af ppMulti. Skal goeres inden kald af ppMulti.

## save_temp: Gem undervejs? Pas p? med parallellisering
## Vigtig ??ndring: Nu skal gr??nser angives korrekt ift. alle parametre, dvs length(upper) = n_par_warp + n_par_amp. De ikke-brugte gr??nser m?? gerne v??re tomme
## Tilsvarende for lower (faktisk nok mere intuitivt). aendringen er kun relevant hvis man ikke maksimerer over alle parametre.



#' Simultaneous inference for Misaligned Multivarariate functional data
#'
#' @param y List of observations in matrix form. NAs allowed
#' @param t List of corresponding time points. NAs not allowed
#' @param basis_fct Basis function for spline
#' @param warp_fct Warp function
#' @param amp_cov Amplitude covariance. Must be on the form function(t, param)
#' @param warp_cov Warp covariance 
#' @param iter two-dimensional integer of maximal number of outer iterations &
#' maximal number of inner iterations per outer iteration.
#' @param w0 Starting values for warp. Should only be used if you have results frmo a previous run.
#' @param amp_cov_par Starting values for amplitude covariance parameters. There are no defaults.
#' @param paramMax Logical vector. Which amplitude parameters to optimise over? Defaults to all parameters.
#' May be overwritten by supplying control parameters.
#' @param like.alt Do not change
#' @param warp_opt If FALSE, warp covariance parameters are kept fixed. 
#' @param like_optim_control List of control options for optimization in outer loop. See details
#' @param use.nlm Use \code{nlm} instead of \code{optim} for optimization? First index for outer loop, second index for inner loop.
#' @param pr Printing option.
#' @param design Design for the experiments. Should be given as a list of one-dimensional vectors or as a design matrix.
#' @param inner_parallel Should the inner optimization be done parallelly?
#' @param save_temp Save estimates after each outer iteration? NULL or the file path.
#' @param functional: Optional and experimnetal parameter that allows the user to use the experimental Functionality package. 
#' Value should be a functional call with one optional argument(ownDeriv = ) that initializes the functional object.
#' If ownDeriv = TRUE the native derivative of the functional object is used, otherwise finite differences are used.
#' 
#'
#' @details ppMulti returns a warning if applied on one-dimensional functional data.
#' 
#' Control parameters:
#' 
#' 
#' lower, upper, method, ndeps, maxit will be sent to optim. See optim for more details. 
#' If use.nlm is selected the optimization is performed using the nlm function. 
#' Bounds are handled through a logit transformation. Note that the optimization will not be able to actually reach the bounds.
#' Presently, the user cannot pass any control values to nlm.
#' 
#' The first entries of lower/upper correspond to warp parameters, while the rest corresponds to
#' amplitude parameters. ppMulti does match upper/lower with corresponding entries in amp_cov, which is important when not all parameters are maximized over. 
#' This is for consistency with randomCycle and optimRule.
#' 
#' randomCycle and optimRule are two ways of optimizing on only a subset of the parameters at a time. These overwrites the paramMax argument.
#' TBD: descriptions of these.
#' 
#' Regarding parallelization: The user must initiate and control parallelization objects by herself. 
#' simm.fda just calls \code{foreach} and does not provide any tools for handling parallelization objects (this is a deliberate design strategy).
#' Note that parallelization does not work in every situation and may depend on operative system among other things.
#'
#'
#' @seealso 
#' Important details can be found in simm-fda-short-desc.pdf
#' 
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
                    w0 = NULL, use.nlm = c(FALSE, FALSE), functional = NULL, 
                    amp_cov_par=NULL, paramMax = rep(TRUE, length(amp_cov_par)),  like.alt = FALSE, warp_opt = TRUE,
                    like_optim_control = list(), pr=FALSE, design = NULL, inner_parallel = FALSE, save_temp = NULL) {

  ## Opsætnings-ting
  
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
    y[[i]] <- y[[i]][!missing_indices, , drop = FALSE] ##Skal være matricer aht. senere
    t[[i]] <- t[[i]][!missing_indices]
  }
 
  t_warped <- t  # Stored warped time
  m <- sapply(y, nrow) # Update m with cleaned data
  
  
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
    require("Functional")
    
    bf0 <- functional()
    basis_fct <- setup.functional.basis(bf0, t[[1]][1], isTRUE(formals(functional)$ownDeriv))
    
    on.exit({ ## Close connection after use.
      print("Closing")
      Functional::removeMember(bf0)
    })
    }
  
  
  
  
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
  if (is.null(w0))   w <- array(attr(warp_fct, 'init'), dim = c(mw, n))
  else w <- w0
  
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
    
    
    #  wis <- list()
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
          gr <- NULL
          w_res <- list()
        
          if (inner_parallel)  w_res <- 
              foreach(i = 1:n) %dopar% {
                if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
                else cis[[i]] <- c 

                warp_optim_method <- 'CG'
                if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate 
                else  ww <- optim(par = w[, i], fn = posterior.lik, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
                
                if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
                return(ww)
              }
          else for ( i in 1:n) {
              
              if (!is.null(design)) cis[[i]] <- c %*%  ( design[[i]] %x% diag(K)  )
              else cis[[i]] <- c 
              warp_optim_method <- 'Nelder-Mead'
              if (use.nlm[2]) ww <- nlm(f = posterior.lik, p = w[,i], warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = cis[[i]], Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$estimate 
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
    
    ## Outer loop part. 
     
    # Construct residual vector for given warp prediction

    Zis <- list()
    r <- y
    
    for (i in 1:n) {
      if (is.null(design)) cis[[i]] <- c
      else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
      
      # Compute warped time
      twarped <- t_warped[[i]] <- warp_fct(w[, i], t[[i]])
      
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        
        Zis[[i]] <- matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, 1]), m[i], mw)
        if (K > 1) for (k in 2:K) {
          Zis[[i]] <- rbind(Zis[[i]], matrix(Zi(twarped, dwarp[[i]], basis_fct, cis[[i]][, k]), m[i], mw))
        }
        
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
        if (warp_opt) {
          param.w <- pars[p_warp]
          par[par1]<- pars[- (p_warp)]
        }
        else {
          param.w <- warp_cov_par
          par[par1]<- pars
        }
        likelihood(par, param.w, r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw, pr = pr)
        
      }
      else like_fct <- function(pars) {
        par <- amp_cov_par
        par[par1]<- pars
        like.S(par, r = r, amp_cov = amp_cov, t = t, pr = pr)
      }
      
      
      # Likelihood gradient
      {
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
        like_optim <- nlm.bound(fct = like_fct , p = paras, lower = lower0, upper = upper0)
        param <- like_optim$estimate
        like_optim$value <- like_optim$minimum
      }
      else {    ## optim optimization
        like_optim <- optim(par = paras, like_fct, gr = like_gr, method = method, lower = lower0, upper = upper0, control = list(ndeps = ndeps, maxit = 20))
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
        
        # Update covariances
        for (i in 1:n) {
          twarped <- t[[i]]
          
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






pavpop_returner_zogr.ny <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5),
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
  
  cat('Inner \t:\tEstimates\n')
    
    wis <- list()
    
    if (ninner > 0) for (iinner in 1:ninner) {
      # Inner loop
      cat(iinner, '\t')
      
      
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
      if (!is.null(design)) {
        c <- (splw.d(y, t, warp_fct, w, Sinv, basis_fct, K = K, design=design))
      } else {
        c <- splw(y, t, warp_fct, w, Sinv, basis_fct, K = K)
      }
      
      
      
      if (warp_change[2] < 3e-3 / sqrt(mw)) break #TODO: Consider other criteria
    }
    
    
    
    
    
    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    
    for (i in 1:n) {
      # compute spline coefficients and warped time
      
      if (is.null(design)) cis[[i]] <- c
      else cis[[i]] <- c %*%  (design[[i]] %x% diag(K))
      twarped <- t_warped[[i]] <- warp_fct(w[, i], t[[i]])
      
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
   
        Zis[[i]] <- multi.Zi(twarped, dwarp[[i]], basis_fct, cis[[i]], mw)
      } else {
        Zis[[i]] <- Matrix(0, m[i]*K, mw)
      }
      
      if (nrow(w) != 1) {
        r[[i]] <- y[[i]] - as.numeric(bf(twarped) %*% cis[[i]]) + as.numeric(Zis[[i]] %*% w[, i])
      }
      
      else {
        rrr <- y[[i]]
        
        for (k in 1:K) {
          rrr[,k] <- rrr[,k] - bf(twarped) %*% cis[[i]][,k] + Zis[[i]][(m[i]*(k-1)+1):(m[i]*k),] * w[,i]
        }
        r[[i]] <- rrr
      }
    }
    
    w_best <- w
    c_best <- c
    
    sigma <- likelihood(amp_cov_par, warp_cov_par, r, Zis, amp_cov, warp_cov, t, tw, sig=T)
    
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best, Zis = Zis, r = r))
}


