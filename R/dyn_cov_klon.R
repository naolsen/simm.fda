

if(F) {

### Eksperimentiel udgave af make_dyn_cov
  
##  Makes use of do.call
  
## All parameters supplied 

make_dyn_cov <- function( timefunction, knots, K, p0) {
  
  indices <- which(lower.tri(matrix(0,K,K), diag = T))
  lpar <- length(indices)
  d <- length(knots)
  cat("This function will have", lpar*d, "parameters \n")
  
  return( function(t, param, ...) {
    a <- list()
    

    
    for (i in 1:d) {
      mat <- matrix(0, K, K)
      mat[indices] <- param[p0 + ((d-1)*lpar+1):(d*lpar) ]
      a[[i]] <- mat %*% t(mat)
    }
    
    b <- alist(t = t, a= a, tw = knots, timefunction = timefunction)
#    b$t <- t
#    b$a <- a
#    b$tw <- knots
    
#    print(timefunction)
    
#    b$timefunction <- timefunction
  #  b[[5]] <- param[1:p0]
    
    print(as.list( param[1:p0]))
    print("fmlksdfmmklsd")
 #   print(formals(timefunction))
    
  #  append(b, as.list( param[1:p0]))
    b <- c(b, as.list( param[1:p0]), ...)
    print(b)
    
    knots
#    kovMat(t, a, knots, timefunction, ...)
    do.call(kovMat, args = b)
    
  })
}


g <- make_dyn_cov(poly.Matern, c(0, 0.4, 1), 4)

g(ti, 1:33, range = 0.1, smooth = 2, koef = c(1,1))




g <- make_dyn_cov(timefunction = poly.Matern, c(0, 0.4, 1), 4, 4)


g(ti, c(range = 0.1, smooth = 2, koef = c(1,1), 1:30))


g <- make_dyn_cov(timefunction = Materntid, c(0, 0.4, 1), 4, 2)


mat <- g(ti, c( range = 0.1,  smooth = 1.2, 1:30), noise = 0)



}
