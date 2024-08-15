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
  
  # Randomly select s of the rows to contain the nonzero rows
  nonzero_indices <- sample(1:p, size=s, replace=FALSE)
  
  
  # Function for sample v times from the set [a,b] U [c,d]
  sample_disjoint <- function(a,b,c,d){
    y <- stats::runif(1, 0, b-a+d-c)
    if( y < (b-a) ){
      x <- a + y
    }else{
      x <- c + y - (b-a)
    }
    return(x)
  }
  
  # Generate the nonzero entries of B0
  for(m in 1:length(nonzero_indices)){
    t <- sample(1:length(response_types), 1)
    nonzero_entries <- sort(sample(1:length(response_types), size=t, replace=F))
      
    if(t==1){
      if(response_types[nonzero_entries]=="continuous" || response_types[nonzero_entries]=="binary") 
        B0[nonzero_indices[m], nonzero_entries] <- sample(c(-1.5,1.5), size=1)
      if(response_types[nonzero_entries]=="count")
        B0[nonzero_indices[m], nonzero_entries] < sample(c(-0.75,0.75), size=1)
    } else if(t>1) {
      
      for(w in 1:length(nonzero_entries)){
        if(response_types[nonzero_entries[w]]=="continuous" || response_types[nonzero_entries[w]]=="binary")
          B0[nonzero_indices[m], nonzero_entries[w]] <- sample_disjoint(-2,-0.5,0.5,2)
        if(response_types[nonzero_entries[w]]=="count")
          B0[nonzero_indices[m], nonzero_entries[w]] <- sample_disjoint(-0.8,-0.4,0.4,0.8)
      }
    }
  }
  return(B0)
}
  
############################
# GENERATE DESIGN MATRIX X #
############################
# n = number of samples
# p = number of covariates
# rho = autocorrelation
# sigma2 = global variance 

generate_X = function(n, p, rho=0.5){
  # Make sure that n and p are valid numbers.
  n <- ceiling(n)
  p <- ceiling(p)
  if(p<=0 || n<=0)
    stop("Error: n and p should both be non-negative.")
  
  # Make sure that rho and sigma2 are suitable values.
  if(rho <= -1 || rho>=1)
    stop("Error: rho should be between -1 and 1 AND sigma2 should be strictly positive.")
  
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  U <- rho^H
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

generate_Y = function(X, B0, response_types, rho=0.5, r=50){
  # Make sure that response_types contains only 'continuous', 'binary', or 'count'
  vals <- unique(response_types)
  if(!all(vals %in% c("binary","count","continuous")))
    stop("Error: response types must be one of: 'binary', 'count', or 'continuous.'")
  
  # Generate random effects u 
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(B0)[2]
  
  # Generate qxq Sigma0 with compound symmetry (CS) structure
  Sigma0 <- matrix(rho, q, q)
  diag(Sigma0) <- rep(1,q)
  mu <- rep(0,q)
  
  # Generate u with signal-to-noise (SNR) of one
  noise <- mvtnorm::rmvnorm(n, mu, Sigma0)
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  Sigma_X <- rho^H
  sigma_e <- sqrt(sum(diag(t(B0)%*%Sigma_X%*%B0))/(sum(diag(t(noise)%*%noise))))
  u <- sigma_e*noise
  
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
