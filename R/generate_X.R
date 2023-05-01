############################
# GENERATE DESIGN MATRIX X #
############################
# n = number of samples
# p = number of covariates
# rho = autocorrelation
# sigma2 = global variance 

generate_X = function(n, p, rho=0.5, sigma2=1){
  # Make sure that n and p are valid numbers.
  n <- ceiling(n)
  p <- ceiling(p)
  if(p<=0 || n<=0)
    stop("Error: n and p should both be non-negative.")
  
  # Make sure that rho and sigma2 are suitable values.
  if(rho <= -1 || rho>=1 || sigma2 <= 0)
    stop("Error: rho should be between -1 and 1 AND sigma2 should be strictly positive.")
  
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  U <- sigma2*rho^H
  mu <- matrix(0, p)
  X <- mvtnorm::rmvnorm(n, mu, U)
  return(X)
}