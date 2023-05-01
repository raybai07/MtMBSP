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