##########################################
# GENERATE TRUE COVARIANCE MATRIX SIGMA0 #
##########################################
# q = number of response variables
# rho = autocorrelation
# sigma2 = global variance 

generate_Sigma0 = function(q, rho=0.5, sigma2=1){
  # Make sure that q is a valid number.
  q <- ceiling(q)
  if(q<=0)
    stop("q should be non-negative.")
  
  # Make sure that rho and sigma2 have suitable values.
  if(rho <= 0 || rho>=1 || sigma2 <= 0)
    stop("Error: rho should be between 0 and 1 AND sigma2 should be strictly positive.")
  
  # Generate qxq Sigma0 with compound symmetry (CS) structure
  Sigma0 <- matrix(rho, q, q)
  diag(Sigma0) <- rep(1,q)
  
  # Return Sigma0
  return(sigma2*Sigma0)
}