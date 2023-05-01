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