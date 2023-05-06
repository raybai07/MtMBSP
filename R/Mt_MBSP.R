#####################################################
# FUNCTION TO IMPLEMENT THE MIXED-TYPE MULTIVARIATE #
# BAYESIAN MODEL WITH SHRINKAGE PRIORS              #
#####################################################

## Authors: Dr. Ray Bai and Dr. Shao-Hsuan Wang
## Email:   raybaistat@gmail.com and picowang@gmail.com 

# INPUT:
# X = n-by-p design matrix
# Y = n-by-q response matrix
# response_types = vector of q response types. The entries must be one of the following: "continuous", "binary", or "count"
# u = first parameter in the TPBN family. Use (u,a)=(0.5,0.5) for the horseshoe prior
# a = second parameter in the TPBN family. Use (u,a)=(0.5,0.5) for the horseshoe prior
# tau = global shrinkage parameter
# d1 = degrees of freedom in inverse-Wishart prior on Sigma
# d2 = scale parameter in inverse-Wishart prior on Sigma
# c1 = shape parameter in gamma prior on dispersion parameter r for count data. Ignored for non-count responses.
# c2 = rate parameter in gamma prior on dispersion parameter r for count data. Ignored for non-count responses.
# algorithm = one-step ("1step") or two-step ("2step")
# step1_iter = number of iterations to run in step 1 of the two-step algorithm. This argument is ignored if
#              algorithm=="1step". If algorithm=="2step", then default is 200  
# bound_error = threshold gamma in Step 1 of the two-step algorithm. 
#               Default of 0 if we run the one-step-algorithm and 0.02 if we run
#               the two-step algorithm.
# max_iter = total number of iterations to run in the Gibbs sampler
# burnin = total number of burn-in iterations for the Gibbs sampler
# details = Boolean variable for whether to return the 0.025 and 0.975 quantiles 
#           and MCMC samples for the model parameters.

# OUTPUT:
# B_est = posterior median estimate for p-by-q regression coefficients matrix B
# B_lower = 0.025 quantile for posteriors of entries in B
# B_upper = 0.975 quantile for posteriors of entries in B
# B_active = binary matrix with "1" for selected variable and "0" for inactive variable
# candidate_set = If two-step method was used, this is the candidate set after Step 1
# B_samples = MCMC samples of B saved after burnin
# Sigma_est = posterior median estimate for Sigma
# Sigma_lower = 0.025 quantile for posterior of Sigma
# Sigma_upper = 0.975 quantile for posterior of Sigma
# Sigma_samples = MCMC samples of Sigma saved after burnin

Mt_MBSP = function(X, Y, response_types, 
                   u=0.5, a=0.5, tau=.001, d1 = dim(Y)[2], d2=10, c1=10, c2=1,
                   algorithm = "1step",
                   step1_iter = 200, bound_error=0, 
                   max_iter = 2000, burnin=1000, 
                   details=TRUE){
  
  # n, p, q
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  
  # Check that X and Y have equal number of samples
  if(dim(Y)[1] != n)
    stop("Error: X and Y must have the same number of samples.")
  # Make sure that response_types vector is of same length as number of responses 
  # and contains only 'continuous', 'binary', or 'count'
  if (q != length(response_types))
      stop("Error: the length of response_types is not equal to the number of responses in Y.")
  vals <- unique(response_types)
  if(!all(vals %in% c("binary","count","continuous")))
    stop("Error: response types must be one of: 'binary', 'continuous', or 'count'.")
  # Check that sample size is sufficiently large.
  if(n < 5*q)
    stop("Error: Estimation is not reasonable if sample size n is too small.")
  # Check that hyperparameters are reasonable.
  if(u <=0 || a <=0)
    stop("Error: TPBN parameters u and a must be strictly positive.")
  if(tau <= 0 || tau>=1)
    stop("Error: Global shrinkage parameter tau should be strictly positive and between 0 and 1.")
  if(d1 <= q-1)
    stop("Error: Inverse-Wishart degrees of freedom should be strictly greater than q-1.")
  if(d2 <= 0)
    stop("Error: Inverse-Wishart scale parameter must be strictly positive.")
  if(c1 <= 0 || c2 <= 0)
    stop("Error: Gamma parameters c1 and c2 must be strictly positive.")
  # Check to ensure that other arguments are reasonable.
  if(max_iter <= burnin || burnin <= 0 || max_iter <= 0)
    stop("Error: Please ensure reasonable arguments for max_iter and burnin, with max_iter > burnin.")
  if(max_iter+burnin<=20)
    stop("Error: max_iter is too small.")
  if(!algorithm %in% c("1step","2step"))
    stop("Error: algorithm must be either '1step' or '2step'.")
  if(algorithm=="2step"){
    if(step1_iter <= 0 || step1_iter >= max_iter || step1_iter >= burnin || step1_iter+burnin >= max_iter)
      stop("Error: If using 2-step algorithm, we require step1_iter < burnin and step1_iter+burnin < max_iter.")
    if(bound_error < 0)
      stop("If using the 2-step algorithm, please specify a bound_error greater than 0.")
    if(is.na(bound_error)){ 
      # Set default bound error of 0.02 for bound_error in the 2-step algorithm if not specified by user
      bound_error = 0.02
    }
  }

  # Initial guesses for Gibbs sampler
  
  # Initial guesses for B and Sigma 
  B_init <- matrix(0, p, q)
  Sigma_init <- diag(q)
  
  # We can refine our initial guesses for B and Sigma for the entries that do NOT
  # correspond to count outcomes.
  which_not_count = which(response_types!='count')
  
  if(length(which_not_count)!=0){
    lambda <- 0.01
    XtY_not_count = t(X) %*% Y[,which_not_count]
    
    # Refine guess for the columns of B that do **not** correspond to count outcomes  
    if(p <= n){
       B_init[,which_not_count] <- chol2inv(chol(t(X)%*%X + lambda*diag(p))) %*% XtY_not_count
    } else if(p > n) { 
       # Using Woodbury matrix identity
       term1 <- (1/lambda)*XtY_not_count
       term2 <- (1/lambda^2)*t(X) %*% chol2inv(chol((1/lambda)*X%*%t(X) + diag(n))) %*% X %*% XtY_not_count
       B_init[,which_not_count] <- term1 - term2
    }
    # Refine guess for the entries of Sigma that do **not** correspond to count outcomes
    resid <- Y[,which_not_count] - X%*%B_init[,which_not_count]
    Sigma_init[which_not_count, which_not_count] <- (n-1)/n * stats::cov(resid)
  }
  
  # Initial guesses for zeta and nu
  zeta_init <- rep(a*tau,p) 
  nu_init <- u*zeta_init 
  
  if(algorithm=="1step"){
    # If running 1-step algorithm, set bound error to be 0 just for good measure
    bound_error <- 0
    # Run Gibbs sampler
    output <- Mt_MBSP_Gibbs(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                            max_iter, burnin, bound_error, details,
                            B_init, Sigma_init, zeta_init, nu_init)
    # Return output
    return(output)
    
  } else if(algorithm=="2step"){
    # If running 2-step algorithm
    
    cat("Stage 1:", "\n")
    # Step 1: Run the Gibbs sampler for step1_iter iterations
    step1_output <- Mt_MBSP_Gibbs(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                                  max_iter=step1_iter, burnin=ceiling(step1_iter/2), 
                                  bound_error, details,
                                  B_init, Sigma_init, zeta_init, nu_init)
    
    # Check the set J of active predictors
    set_J <- which(rowSums(step1_output$B_active) != 0)
    
    # If Stage 1 result 
    if(length(set_J)==0){
         cat("Stage 1 resulted in a null model. Consider increasing bound_error.", "\n")
         return(step1_ouput)
    } else if(length(set_J) >= n) {
        # If the set J has n or more predictors, further reduce its size to n-1
        q_j <- rep(0, length(set_J))
        for(jj in 1:length(q_j)){
          q_j[jj] <- max(abs(step1_output$B_est[jj, ])) 
        }
        max_q_j_indices <- sort(utils::tail(order(q_j), n-1))
        set_J <- set_J[max_q_j_indices]
    }
    
    # Estimated set of inactive predictors  
    set_Jc <- setdiff(seq(1:p), set_J)

    cat("\n", "Stage 2:", "\n")
    # Step 2: Run the Gibbs sampler with the smaller subset J
    # Estimated set of active predictors
    
    # Initialize B, Sigma, zeta_init, and nu_init for Gibbs sampler  
    X_2 <- as.matrix(X[, set_J]) # Step 2 is run with only the variables in candidate set J_n
    B_init_2 <- step1_output$B_est[set_J,]
    if(dim(X_2)[2]==1){ B_init_2 <- t(B_init_2) }
    Sigma_init_2 <- step1_output$Sigma_est
    zeta_init_2 <- rep(a*tau,length(set_J)) 
    nu_init_2 <- u*zeta_init 
    
    # New max_iter and burnin
    step2_iter = max_iter - step1_iter
    step2_burn = burnin - step1_iter 
    
    # [1] Run Step 2 of the Gibbs sampler for max_iter-step1_iter iterations
    # and no threshold gamma (bound_error=0)
    step2_output <- Mt_MBSP_Gibbs(X=X_2, Y, response_types, u, a, tau, d1, d2, c2, c2,
                                  max_iter=step2_iter, burnin=step2_burn, 
                                  bound_error=0, details=details,
                                  B_init=B_init_2, Sigma_init=Sigma_init_2, 
                                  zeta_init=zeta_init_2, nu_init=nu_init_2)
    # New output
    B_est <- matrix(0, p, q)
    B_active <- matrix(0, p, q)
    B_lower <- matrix(0, p, q)
    B_upper <- matrix(0, p, q)
    
    B_est[set_J, ] <- step2_output$B_est
    B_active[set_J, ] <- step2_output$B_active
    B_lower[set_J, ] <- step2_output$B_lower
    B_upper[set_J, ] <- step2_output$B_upper
    B_lower[set_Jc, ] <- -1e-5
    B_upper[set_Jc, ] <- 1e-5
      
    if(details==TRUE){
      output <- list(
        B_est = B_est,
        B_active = B_active,
        B_lower = B_lower,
        B_upper = B_upper,
        candidate_set = set_J,  
        B_samples = step2_output$B_samples,
        Sigma_est = step2_output$Sigma_est,
        Sigma_lower = step2_output$Sigma_lower,
        Sigma_upper = step2_output$Sigma_upper,
        Sigma_samples = step2_output$Sigma_samples
      )
    } else{
      list(
        B_est = B_est,
        B_active = B_active,
        Sigma_est = step2_output$Sigma_est)    
    }
  }
}


##########################
# GIBBS SAMPLER FUNCTION #
##########################

# Gibbs sampler for Mt-MBSP
Mt_MBSP_Gibbs = function(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                         max_iter, burnin, bound_error, details,
                         B_init, Sigma_init, zeta_init, nu_init){

  # sigmoid function
  sigmoid <- function(x){(1+exp(-x))^{-1}}

  # Extract dimensions n, p, and q
  # X is a n times p times n covariate matrix
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  
  ################## 
  # Initialization #
  ##################
 
  # Initialize Gibbs sampler
  B <- B_init
  Sigma <- Sigma_init
  zeta <- zeta_init
  nu <- nu_init
  # Initial guess for latent matrices Z and W
  Z <- Y
  W <- matrix(1,n,q) # default for Gaussian
  # Initial guess for latent matrix latentU ~ N(0, Sb);
  latentU <- Y-X%*%B
  
  # Matrix for updating the dispersion parameter r for count responses.
  # To make it consistent with response_types, only the indices that correspond
  # to count variables are updated, the rest are NA.
  dispersion_r <- rep(NA, q)
  # Initial guess for dispersion parameters
  dispersion_r[which(response_types=="count")] <- 10
    
  #################
  # Gibbs sampler #
  #################
  # Lists to hold the draws of B of Sigma
  B_samples <- rep(list(matrix(0,p,q)), max_iter)
  Sigma_samples <- rep(list(matrix(0,q,q)), max_iter)
  
  for(bj in 1:max_iter){
    if (bj %% 100 == 0)
      	cat("Iteration:", bj,"/",max_iter, "\n")
    
    ## Update nxq polya-Gamma weight matrix W

    for(j in 1:q){
      thetaj <- (X%*%B[,j]+latentU[,j])
      W[,j] <- switch(response_types[j],
                      binary = BayesLogit::rpg(n,1, thetaj),
                      continuous = rep(1,n),
                      count = BayesLogit::rpg(n,Y[,j]+dispersion_r[j], thetaj))
    
      # For count responses, update dispersion parameters
      
      if(response_types[j]=='count'){
        # Update dispersion parameter   
        pp <- as.vector(sigmoid(thetaj))
        pp[pp>.9999] <- .9999
        
        # Initialize l vector
        l <- rep(0,n)
        # Fill in the l vector
        for(ii in 1:n){
          pj <- round(dispersion_r[j]/(dispersion_r[j]+1:Y[ii,j]-1), 6)
          l[ii] <- sum(stats::rbinom(Y[ii,j], 1, pj))
        } 
        # Update r from conjugate gamma distribution given l and beta
        dispersion_r[j] <- stats::rgamma(1, c1+sum(l), c2-sum(log(1-pp)))
      }
    }
    vecW <- as.vector(W) #reset W.   
  
    # Update latent matrix Z
    # for Z = (Y-U)/W; W==1 if Y is continuous
  
    for(j in 1:q){
      Z[,j] <- switch(response_types[j],
                      binary = (Y[,j]-1/2)/W[,j],
                      continuous = Y[,j],
                      count = (Y[,j]-(Y[,j]+dispersion_r[j])/2)/W[,j])
    }

    # Update B
    for(j in 1:q){
      if(p<=n){
         if(length(zeta)>1){
           Delta_inv <- chol2inv(chol(t(X)%*%(X*W[,j])+diag(1/zeta)))
         } else if(length(zeta)==1){
           Delta_inv <- chol2inv(chol(t(X)%*%(X*W[,j])+1/zeta))
         }
         post_M <- Delta_inv %*% t(X*W[,j]) %*%(Z[,j]-latentU[,j])
         B[,j] <- mvtnorm::rmvnorm(1, post_M, Delta_inv)
      }
      else if(p>n){
        # We use the fast sampling method from Bhattacharya et al. (2016)
        post_M <- t(matrix(W[,j]*(Z[,j]-latentU[,j]),1,n)%*%X)
        zU <- matrix(mvtnorm::rmvnorm(1,rep(0,n),diag(W[,j])),n,1);
        zM <- stats::rnorm(p,0,1)
        
        # formula diag(x)%*%H=x*H for px1 vector x and p x n matrix H. 
        Jm <- zeta*(post_M-(t(X)%*%zU+(1/sqrt(zeta))*zM))
        JW <- Jm- (zeta*t(X))%*%(chol2inv(chol(diag(1/W[,j])+X%*%(zeta*t(X))))%*% (X %*%Jm))
        B[,j] <- JW
      }
    } 

    # Update Sigma
    Sigma <- MCMCpack::riwish(n+d1, (t(latentU)%*%latentU)+diag(q)*d2)
  
    # Update latent variables matrix U
    tmp_mat <- Z-X%*%B
    for(i in 1:n){
       Omega_i <- diag(W[i,])
       mm <- Omega_i%*%tmp_mat[i,]
       JJ <- chol2inv(chol(Omega_i+Sigma))
       latentU[i,] <- mvtnorm::rmvnorm(1, JJ%*%mm, JJ)
    }

    # Update zeta_i's and nu_i's
    for (i in 1:p){
      norm_term <- sum((t(B[i,]))^2) 
      v <- max(norm_term, .Machine$double.eps) # to prevent chi parameter from collapsing to 0
      zeta[i] <- GIGrvg::rgig(n=1, lambda=u-q/2, chi=v, psi=2*nu[i])
      nu[i] <- stats::rgamma(n=1, shape=a, scale=1/(tau+zeta[i]))
    }
    
    B_samples[[bj]] <- B;
    Sigma_samples[[bj]] <- Sigma;
  }
  
  ########################
  # end of Gibbs sampler #
  ########################
  
  # Discard burn-in 
  B_samples <- utils::tail(B_samples, max_iter-burnin) # only contain the tail part
  Sigma_samples <- utils::tail(Sigma_samples, max_iter-burnin) # only contain the tail part

  # Extract the posterior median and 2.5th and 97.5th quantiles of B 
  arrB<- array(unlist(B_samples), c(p,q,length(B_samples)))
  B_est<-NULL; B_lower<-NULL;B_upper<-NULL;active_predictors<-NULL
  for(j in 1:q){
    mbsp_quantiles <- apply(arrB[,j,], 1, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
    
    # Take posterior median as point estimate for B
    B_est <- cbind(B_est, mbsp_quantiles[2,])
    # For marginal credible intervals
    B_lower <- cbind(B_lower, mbsp_quantiles[1,])
    B_upper <- cbind(B_upper, mbsp_quantiles[3,])
    
  }
  
  # Set J
  B_lower <- ifelse(B_lower > -bound_error*stats::sd(B_lower), bound_error*stats::sd(B_lower), B_lower)
  B_upper <- ifelse(B_upper < bound_error*stats::sd(B_upper), -bound_error*stats::sd(B_upper), B_upper)
  # for active.predictors 
  temp <- sign(B_lower)*sign(B_upper)
  active_predictors <- 1*(temp>=0)

  # Return list of B_est, B_lower, 
  # and B_upper, active.predictors 

  B_est <- matrix(B_est, p, q, byrow=F)
  B_lower <- matrix(B_lower, p, q, byrow=F)
  B_upper <- matrix(B_upper, p, q, byrow=F)
  
  # Extract the posterior median and 2.5th and 97.5th quantiles of Sigma 
  arrS<- array(unlist(Sigma_samples), c(q,q,length(Sigma_samples)))
  Sigma_est<-NULL; Sigma_lower<-NULL;Sigma_upper<-NULL;
  for(j in 1:q){
    Smbsp_quantiles <- apply(arrS[,j,], 1, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
    
    # Take posterior median as point estimate for B
    Sigma_est <- cbind(Sigma_est, Smbsp_quantiles[2,])
    # For marginal credible intervals
    Sigma_lower <- cbind(Sigma_lower, Smbsp_quantiles[1,])
    Sigma_upper <- cbind(Sigma_upper, Smbsp_quantiles[3,])
  }
  
  # Return list of Sigma_est, Sigma_lower, Sigma_upper
  Sigma_est <- matrix(Sigma_est, q, q, byrow=F)
  Sigma_lower <- matrix(Sigma_lower, q, q, byrow=F)
  Sigma_upper <- matrix(Sigma_upper, q, q, byrow=F)
  

  if(details==TRUE){
    output <- list(B_est = B_est,
                   B_active = active_predictors,
                   B_lower = B_lower,
                   B_upper =  B_upper,
                   B_samples = B_samples,
                   Sigma_est = Sigma_est,
                   Sigma_lower = Sigma_lower,
                   Sigma_upper = Sigma_upper,
                   Sigma_samples = Sigma_samples)
  } else{
    output <- list(B_est = B_est,
                   B_active = active_predictors,
                   Sigma_est = Sigma_est)
  }
  
  # Return list
  return(output)
}
