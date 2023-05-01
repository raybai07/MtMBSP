###########################################
# FUNCTIONS FOR GENERATING SYNTHETIC DATA #
###########################################
# Scripts to generate synthetic datasets from Section 5.1 of the paper

## Authors: Dr. Ray Bai and Dr. Shao-Hsuan Wang
## Email:   raybaistat@gmail.com and picowang@gmail.com 


############################
# GENERATE TRUE REGRESSION #
# COEFFICIENTS MATRIX B0   #
############################
# p = number of rows in B0
# s = number of nonzero rows
# response_types = vector of q response types ("continuous", "binary", and "count")

generate_B0 = function(p, s, response_types){
  # Make sure that p and s are valid numbers
  p <- ceiling(p)
  s <- ceiling(s)
  if(p<=0 || s<=0 || s>p)
    stop("Error: p and s should be positive and p<=s.")
  
  # Make sure that response_types contains only 'continuous', 'binary', or 'count'
  vals <- unique(response_types)
  if(!all(vals %in% c("binary","count","continuous")))
    stop("Error: response types must be one of: 'binary', 'continuous', or 'count'.")
  
  # Initialize matrix B0
  q <- length(response_types) 
  B0 <- matrix(0,p,q)
  
  # Randomly select 10 of the row indices to contain the nonzero rows
  nonzero_indices <- sample(1:p, size=s, replace=FALSE)
  
  # Function for sample n times from the set [a,b] U [c,d]
  sample_disjoint <- function(n, a,b,c,d){
    # Initialize
    sample_vec <- rep(0,n)
    
    # Sample entries in sample_vec
    for(ntimes in 1:n){
      y <- stats::runif(1, 0, b-a+d-c)
      if( y < (b-a) ){
        x <- a + y
      }else{
        x <- c + y - (b-a)
      }
      sample_vec[ntimes] <- x
    }
    return(sample_vec)
  }
  
  # Generate the nonzero entries of B0
  for(k in 1:length(response_types)){
    if(response_types[k]=="continuous" || response_types[k]=="binary"){
      B0[nonzero_indices,k] <- sample_disjoint(s,-5,-0.5,0.5,5)
    } else if(response_types[k]=="count"){
      # Sample either -1 or 1
      sgn <- sample(c(-1,1), size=1, replace=FALSE)
      if(sgn == 1){
        B0[nonzero_indices,k] <- stats::runif(s, 0.3, 0.6)
      } else if(sgn == -1){
        B0[nonzero_indices,k] <- stats::runif(s, -0.6, -0.3)
      }
    }
  }
  return(B0)
}
  
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

##############################
# GENERATE RESPONSE MATRIX Y #
##############################
# X = design matrix
# B0 = true regression coefficients matrix
# Sigma0 = true covariance matrix
# response_types = vector of q response types ("continuous", "binary", and "count")
# r = number of successful trials. Only for count responses. 
#     Default is r=50 for equidispersed data. Decrease r for overdispersed data. 

generate_Y = function(X, B0, Sigma0, response_types, r=50){
  # Make sure that response_types contains only 'continuous', 'binary', or 'count'
  vals <- unique(response_types)
  if(!all(vals %in% c("binary","count","continuous")))
    stop("Error: response types must be one of: 'binary', 'count', or 'continuous.'")
  
  # Make sure that the matrices are conformable
  q <- length(response_types)
  if(dim(B0)[2]!= q || dim(Sigma0)[1]!=q)
    stop("Error: B0 and Sigma0 must have the same number of columns as the length of response_types.")
  if(dim(X)[2] != dim(B0)[1])
    stop("Error: The number of columns in X should be equal to the number of rows in B0.")
  
  # Generate random effects u
  n <- dim(X)[1]
  u <- mvtnorm::rmvnorm(n,rep(0,q), Sigma0) 

  # Generate theta
  theta <- X %*% B0 + u
  
  # Sigmoid function
  sigmoid <-function(x){1/(1+exp(-x))}  
  
  # Initialize Y
  Y <- matrix(0,n,q)
  # Fill in the columns of Y
  for(k in 1:length(response_types)){
    if(response_types[k]=="continuous"){
       Y[,k] <- matrix(X%*%B0[,k],n,1)+ u[,k]+stats::rnorm(n,0,1)       
    } else if(response_types[k]=="binary"){
        prob <- sigmoid(matrix(X%*%B0[,k],n,1)+ u[,k])
        Y[,k] <- stats::rbinom(n, 1, prob)
    } else if(response_types[k]=="count"){
        prob <- sigmoid(matrix(X%*%B0[,k],n,1)+ u[,k])
        Y[,k] <- stats::rnbinom(n, size = r, prob = 1-prob)
    }
  }
  return(Y)
}