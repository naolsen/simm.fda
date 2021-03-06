

if (F) {
#### Example.R
# Data taken from <MOCAP> ets.

## Data is provided as 
# MCD.data y-values
# MCD.time time points

# The time points have been scaled to be within [0.1, 0.9] ...

# Make basis function
# bf <- make_basis_fct_old(seq(0, 1, len=32)[2:31], intercept=T, boundary = c(0, 1))
bf <- make_basis_fct(seq(0, 1, len=32), intercept=T, type = 'B-spline')

## Make warp function
tw <- seq(0, 1, length = 5)[2:4] ## anchor points for hyman spline
wf.noshift <- make_warp_fct(type="smooth", tw=tw)

## warp function with shift
wf <- w.shift(wf.noshift, 0.25)

warp_cov_par <- c(tau = 0.5)
#warp_cov <- make_cov_fct(Brownian, noise = FALSE, param = warp_cov_par, type = 'bridge')

warp_cov <- default_warp_cov(0.5)
warp_cov <- warp_and_shift_cov(c(0.5, 1))


## For this example we assume no amplitude cross-correlation, otherwise ...
# Matern covariance; smooth parameter 2, unknown range.

wrapper <- function(t, par) {
  simm.fda:::mvMatern(t, par[1], 2, rep(1,3), diag(par[2:4]))
}

lower <- c(1e-4, 1e-4,rep(1e-4, 4))
upper <- rep(1e5, 6)

mcd.res <-ppMulti(MCD.data, MCD.time, bf, wf, amp_cov = wrapper, warp_cov, amp_cov_par = c(0.2, rep(100, 3)), pr = T, paramMax = rep(T, 4),
                  like_optim_control = list(lower = lower, upper = upper), iter = c(31, 7))


}

