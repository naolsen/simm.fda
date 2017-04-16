


### Functionality
# NOT to be exported

setup.functional.basis <- function(addr, some_value, ownDeriv) {
  ## NOT to be used from outside the package
  
  if (ownDeriv) {
    cat('functional form with own derivative \n')
    basis_fct <- function(t, deriv = FALSE) {
      if (deriv) eval.deriv.direct(addr, t, checkValidity = FALSE)
      else eval.direct(addr, t, checkValidity = FALSE)
    }
  }
  else {
    basis_fct <- function(t, deriv = FALSE) { ## Denne konstruktion er mest for konsistens skyld. Bør forbedres
      
      if (deriv) ## Vedtagelse: epsilon = 1e-5
        (eval.direct(addr, t + 1e-5, checkValidity = FALSE) -
           eval.direct(addr, t - 1e-5, checkValidity = FALSE)) / (2e-5)
      
      else eval.direct(addr, t, checkValidity = FALSE)
    }
  }
  
  ## Find size of basis of function (should be changed in the future)
  attr(basis_fct, 'df') <- length(basis_fct(some_value ))
  
  basis_fct
  
}
