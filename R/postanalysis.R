





#' Evaluated warp functions
#'
#' @description Returns and optionally plots predicted warps.
#' 
#' @param w warp values
#' @param warp_fct
#' @param tout Evaluation points. If missing, equidistant evaluation on [0,1] according to seq.length is used in the evaluation
#' @param seq.length number of time points used in the evaluation.
#' @param plot Plot the result? Uses matplot.
#' @param ... arguments to be passed to matplot
#'
#' @return A matrix of dim seq.length x n
#' @export
#'
#' @examples  
#' \donttest{
#' 
#' tw <- c(0.25, 0.5, 0.75)
#' wf <- <- make_warp_fct(type="smooth", tw=tw)
#' 
#' <run estimation on data set and return object res>
#' 
#' predict_warps(res$w, wf, plot = T, type="l", main = "Predicted warp functions")
#' 
#' }
#' 

predict_warps <- function(w, warp_fct, seq.length = 101, plot = TRUE, tout , ...) {
  
  mw <- attr(warp_fct, 'mw')
  if (nrow(w) != mw) stop("w and warp function doesn't match")
  
  n <- ncol(w)
  
  if (!missing(tout)) ti <- tout
  else ti <- seq(0, 1, length = seq.length)
  
  warp.val <- matrix(nc = n, nr =length(ti))

  
  for (i in 1:n) {
    warp.val[,i] <- warp_fct(w[,i], ti)
  }
  if (plot) matplot(ti, warp.val, ...)
  
  return(warp.val)
}




#' Plot aligned data
#'
#' @param y,t data
#' @param w warp. If missing, no alignment is done.
#' @param warp_fct 
#' @param coords Which coordinates of trajectories? E.g. 1 for fist column.
#' @param type plot type
#' @param ... Further arguments to be passed to plot and points. 
#' 
#' @description Much improvement needs to be done on this function.
#' 
#'
#' @return
#' @export
#'
#' 
align_data_plot <- function(y, t, w, warp_fct, coords, type="l", ...) {
  ## At blive gjort.
  2+2

  tw <- attr(warp_fct, 'tw')  
  n <- length(y)
  if (missing(w)) w <- matrix(0, length(tw), n)
  
  if (ncol(w) != n || length(t) != n) stop("Lengths do not match!")

  
  yz <- sapply(y, function (x) x[,coords])
  range <- range(unlist(yz), na.rm = T)
  plot(c(0,1), range, type = "n", ...)
  
  for (i in 1:n) {
    
    wt <- warp_fct(w[,i], t[[i]])
    for (j in coords) {
      points(wt, y[[i]][,j] , type = type, ...)
    }
  }
  
}



