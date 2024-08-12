#####################################################
# FUNCTION TO IMPLEMENT THE MIXED-TYPE MULTIVARIATE #
# BAYESIAN MODEL WITH SHRINKAGE PRIORS              #
#####################################################

library(mvtnorm)
library(GIGrvg) 
library(MCMCpack) 
library(BayesLogit) 
library(stats) 
library(utils)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)

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
# algorithm = one-step ("1step") or two-step ("2step"). Default is "1-step."
# niter = number of MCMC iterations to run for one-step algorithm.
# burn = number of MCMC iterations to discard as burn-in for one-step algorithm.
# step2_niter = number of MCMC iterations to run in step 2 of the two-step algorithm.
# step2_burn = number of MCMC iterations to discard as burn-in in step 2 of the two-step algorithm.
# threshold = threshold gamma to be used in Step 1 of the two-step algorithm. 
#             If the user specifies a single value, then this threshold is used for gamma.
#             If the user specifies a grid, then the function searches for the optimal gamma that minimizes WAIC.
#             Default is to search a grid from 0.02 to 0.40.
#             This argument is ignored if algorithm="1step"
# parallelize = Boolean variable for whether or not to parallelize Step 2 of the two-step algorithm.
#               This argument is ignored if algorithm="1step"
# ncores = number of cores to use for parallelization if parallelize=TRUE

# OUTPUT:
# B_est = posterior median estimate for p-by-q regression coefficients matrix B
# B_active = binary matrix with "1" for selected variable and "0" for inactive variable
# B_lower = lower endpoints of 95% uncertainty intervals for entries in B. 
#           If the two-step algorithm is used and the estimated entry of B_est is exactly zero, then 
#           this returns the lower endpoint of the posterior credible interval obtained from Step 1.
# B_upper = upper endpoints of 95% uncertainty intervals for entries in B.
#           If the two-step algorithm is used and the estimated entry of B_est is exactly zero, then 
#           this returns the upper endpoint of the posterior credible interval obtained from Step 2.
# Sigma_est = posterior median estimate for Sigma
# Sigma_lower = lower endpoints of 95% uncertainty intervals for entries in Sigma
# Sigma_upper = upper endpoints of 95% uncertainty intervals for entries in Sigma
# opt_threshold = optimal gamma chosen by minimizing WAIC. Only returned if "2step" is used for algorithm
# set_J = the initial set J_n chosen in Step 1 if "2step" is used for algorithm
# B_samples = MCMC samples of B saved after burn.
#             If algorithm="2step", then the MCMC samples from Step 1 and Step 2 are both returned for B
# Sigma_samples = MCMC samples of Sigma saved after burn


###################
## MAIN FUNCTION ##
###################
Mt_MBSP = function(X, Y, response_types, 
                   u=0.5, a=0.5, tau=1/(dim(X)[2]*sqrt(dim(X)[1]*log(dim(X)[1]))), 
                   d1 = dim(Y)[2], d2=10, c1=10, c2=1,
                   algorithm = "1step",
                   niter = 1100, burn=100,
                   step2_niter = 1100, step2_burn=100, 
                   threshold = seq(from=0.02, to=0.40, by=0.02),
                   parallelize = TRUE, ncores = 10){
  
  ####################
  ## Error handling ##
  ####################
  
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
  if(tau < 1e-5) tau <- 1e-5 # To prevent tau from being too small
  if(d1 <= q-1)
    stop("Error: Inverse-Wishart degrees of freedom should be strictly greater than q-1.")
  if(d2 <= 0)
    stop("Error: Inverse-Wishart scale parameter must be strictly positive.")
  if(c1 <= 0 || c2 <= 0)
    stop("Error: Gamma parameters c1 and c2 must be strictly positive.")
  # Check to ensure that other arguments are reasonable.
  if(niter <= burn || burn <= 0 || niter <= 0)
    stop("Error: Please ensure reasonable arguments for niter and burn, with niter > burn.")
  if(!algorithm %in% c("1step","2step"))
    stop("Error: algorithm must be either '1step' or '2step'.")
  if(algorithm=="2step"){
    if(step2_niter <= step2_burn || step2_burn <= 0 || step2_niter <= 0)
      stop("Error: Please ensure reasonable arguments for step2_niter and step2_burn, 
            with step2_niter > step2_burn.")
    if(!all(threshold>=0))
      stop("If using the 2-step algorithm, the threshold must be nonnegative.")
    if(all(is.na(threshold))){ 
      # Set default as an equispaced grid from 0.02 to 0.4
      threshold = seq(from=0.02, to=0.4, by=0.02)
    }
    if(length(threshold)==1)
      parallelize = FALSE
    if(ncores <= 0)
      stop("Number of cores should be greater than 0.")
    if((ncores > parallel::detectCores()-1)) 
      ncores <- parallel::detectCores()-1
    ncores <- as.integer(ncores)
  }
  
  #######################
  ## Fit Mt-MBSP model ##
  #######################
  
  if(algorithm=="1step"){
    # Run Gibbs sampler
    output <- Mt_MBSP_Gibbs(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                            niter, burn, nugget=0.03, return_WAIC = FALSE)
  
    # Extract summary statistics
    posterior_summaries <- Mt_MBSP_summary(output$B_samples, output$Sigma_samples,
                                           threshold = 0)
    # Return results
    return(list(B_est = posterior_summaries$B_est,
                B_active = posterior_summaries$B_active,
                B_lower = posterior_summaries$B_lower,
                B_upper = posterior_summaries$B_upper,
                Sigma_est = posterior_summaries$Sigma_est,
                Sigma_lower = posterior_summaries$Sigma_lower,
                Sigma_upper = posterior_summaries$Sigma_upper,
                B_samples = output$B_samples,
                Sigma_samples = output$Sigma_samples))
    
  } else if(algorithm=="2step"){
    
    cat("Running Stage 1.", "\n")
    
    # Step 1: Run the Gibbs sampler for niter iterations
    step1_output <- Mt_MBSP_Gibbs(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                                  niter=niter, burn=burn, nugget=0.02, return_WAIC = FALSE)
    
    # Number of samples to save
    nsave = niter-burn
    step1_B_samples <- utils::tail(step1_output$B_samples, nsave)
    step1_Sigma_samples <- utils::tail(step1_output$Sigma_samples, nsave)
    step1_summaries <- Mt_MBSP_summary(step1_B_samples, step1_Sigma_samples,
                                       threshold = 0)
    
    # To hold results
    threshold <- sort(threshold)
    set_J <- vector(mode = "list", length = length(threshold))
    
    for(l in 1:length(threshold)){
      # Get the set J_n of active predictors
      tmp_summaries <- Mt_MBSP_summary(step1_B_samples, step1_Sigma_samples,
                                       threshold=threshold[l])
      set_J[[l]] <- which(rowSums(tmp_summaries$B_active) != 0)
      
      if(length(set_J[[l]]) >= n) {
        # If the active set has n or more predictors, further reduce its size to n-1
        q_j <- rep(0, length(set_J[[l]]))
        for(jj in 1:length(q_j)){
          q_j[jj] <- max(abs(tmp_summaries$B_est[jj, ])) 
        }
        max_q_j_indices <- sort(utils::tail(order(q_j), n-1))
        set_J[[l]] <- set_J[[l]][max_q_j_indices]
      }
    }
    
    if(length(unlist(set_J))==0){
        cat("Stage 1 resulted in a null model. Consider increasing the minimum threshold.", "\n")
        tmp_summaries <- Mt_MBSP_summary(step1_B_samples, step1_Sigma_samples,
                                         threshold=threshold[1])
        return(list(B_est = matrix(0,p,q),
                    B_active = matrix(0,p,q),
                    B_lower = tmp_summaries$B_lower,
                    B_upper = tmp_summaries$B_upper,
                    Sigma_est = tmp_summaries$Sigma_est,
                    Sigma_lower = tmp_summaries$Sigma_lower,
                    Sigma_upper = tmp_summaries$Sigma_upper,
                    B_samples = step1_B_samples,
                    Sigma_samples = step1_Sigma_samples))
    } 
    
    # Remove the empty elements of set_J 
    nonempty_sets <- as.numeric(lapply(set_J, length) != 0)
    threshold_nonempty <- threshold[which(nonempty_sets!=0)]
    set_J <- set_J[which(nonempty_sets!=0)]
    
    # Run Stage 2 
    if(parallelize==TRUE){
      cat("Running Stage 2 in parallel for each threshold.", "\n")
    
      clusters <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(clusters)
      doRNG::registerDoRNG(1)
      
      step2_samples <- foreach::foreach(set_J=set_J, 
                                        .export=c("Mt_MBSP_summary", "Mt_MBSP_Gibbs"),
                                        .combine=list,
                                        .multicombine=TRUE) %dorng% {
      
        # Step 2: Run the Gibbs sampler with the smaller subset J
        X_2 <- as.matrix(X[, set_J]) # Step 2 is run with only the variables in candidate set J_n
      
        # Estimated set of active predictors
        step2_gibbs <- Mt_MBSP_Gibbs(X=X_2, Y=Y, response_types=response_types, 
                                     u=u, a=a, tau=tau, d1=d1, d2=d2, c1=c1, c2=c2,
                                     niter=step2_niter, burn=step2_burn,
                                     nugget=0.03, return_WAIC=TRUE)
      
        list(B_samples = step2_gibbs$B_samples,
             Sigma_samples = step2_gibbs$Sigma_samples,
             WAIC = step2_gibbs$WAIC,
             set_J = set_J)
      }
      parallel::stopCluster(clusters)
    
    } else {
      
      # Do not parallelize
      step2_samples <- vector(mode = "list", length = length(set_J))
      
      # Step 2: Run the Gibbs sampler with the smaller subset J
      for(k in 1:length(step2_samples)){
        
        cat("Running Stage 2 for threshold", threshold_nonempty[k], "\n")
        
        
        X_2 <- as.matrix(X[, set_J[[k]]]) # Step 2 is run with only the variables in candidate set J_n
        
        # Estimated set of active predictors
        step2_gibbs <- Mt_MBSP_Gibbs(X=X_2, Y=Y, response_types=response_types, 
                                     u=u, a=a, tau=tau, d1=d1, d2=d2, c1=c1, c2=c2,
                                     niter=step2_niter, burn=step2_burn,
                                     nugget=0.03, return_WAIC=TRUE)
        
        step2_samples[[k]] <- list(B_samples = step2_gibbs$B_samples,
                                   Sigma_samples = step2_gibbs$Sigma_samples,
                                   WAIC = step2_gibbs$WAIC,
                                   set_J = set_J)
      }
    }

    # Extract posterior summaries
    nonempty_card <- length(step2_samples)
    step2_summaries <- vector(mode = "list", length = nonempty_card)
    WAIC_vec <- rep(0, nonempty_card)
    
    for(k in 1:nonempty_card){
      step2_summaries[[k]] <- Mt_MBSP_summary(step2_samples[[k]]$B_samples, 
                                              step2_samples[[k]]$Sigma_samples,
                                              threshold=0)
    
      # Store WAIC with threshold[k]
      WAIC_vec[k] <- step2_samples[[k]]$WAIC
    }
    
    # Final model chosen from WAIC
    opt_index <- which.min(WAIC_vec)
    opt_threshold <- threshold_nonempty[opt_index]
    step2_final_samples <- step2_samples[[opt_index]]
    step2_final_summaries <- step2_summaries[[opt_index]]
    set_J_final <- set_J[[opt_index]]
    set_Jc_final <- setdiff(seq(1:p), set_J_final)
    
    # New summaries
    step2_B_est <- matrix(0, p, q)
    step2_B_active <- matrix(0, p, q)
    step2_B_lower <- matrix(0, p, q)
    step2_B_upper <- matrix(0, p, q)
    
    # Estimates
    step2_B_est[set_J_final, ] <- step2_final_summaries$B_est
    step2_B_active[set_J_final, ] <- step2_final_summaries$B_active
    step2_B_lower[set_J_final, ] <- step2_final_summaries$B_lower
    step2_B_lower[set_Jc_final, ] <- step1_summaries$B_lower[set_Jc_final, ]
    step2_B_upper[set_J_final, ] <- step2_final_summaries$B_upper
    step2_B_upper[set_Jc_final, ] <- step1_summaries$B_upper[set_Jc_final, ]
    step2_Sigma_est <- step2_final_summaries$Sigma_est
    step2_Sigma_lower <- step2_final_summaries$Sigma_lower
    step2_Sigma_upper <- step2_final_summaries$Sigma_upper
    
    return(list(B_est = step2_B_est,
                B_active = step2_B_active,
                B_lower = step2_B_lower,
                B_upper = step2_B_upper,
                Sigma_est = step2_Sigma_est,
                Sigma_lower = step2_Sigma_lower,
                Sigma_upper = step2_Sigma_upper,
                opt_threshold = opt_threshold,
                set_J = set_J_final,
                step1_B_samples = step1_B_samples,
                step2_B_samples = step2_final_samples$B_samples,
                Sigma_samples = step2_final_samples$Sigma_samples))
  }
}


############################
## GIBBS SAMPLER FUNCTION ##
############################

# Gibbs sampler for Mt-MBSP
Mt_MBSP_Gibbs = function(X, Y, response_types, u, a, tau, d1, d2, c1, c2,
                         niter, burn, nugget, return_WAIC){

  # sigmoid function
  sigmoid <- function(x){(1+exp(-x))^{-1}}

  # Extract dimensions n, p, and q
  # X is a n times p times n covariate matrix
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  #################### 
  ## Initialization ##
  ####################
 
  # Initialize Gibbs sampler
  B <- matrix(0,p,q)
  Sigma <- diag(q)
  zeta <- rep(1,p)
  nu <- rep(1,p)
  # Initialize W and Z
  W <- matrix(0,n,q)  
  Z <- matrix(0,n,q)
  # Initial guess for latent matrix latentU ~ N(0, Sb);
  latentU <- matrix(0,n,q)
  
  # Matrix for updating the dispersion parameter r for count responses.
  # To make it consistent with response_types, only the indices that correspond
  # to count variables are updated, the rest are NA.
  dispersion_r <- rep(NA, q)
  # Initial guess for dispersion parameters
  dispersion_r[which(response_types=="count")] <- 10
    
  ###################
  ## Gibbs sampler ##
  ###################
  
  # Lists to hold the draws of B, Sigma, W, U, and Z
  B_samples <- rep(list(matrix(0,p,q)), niter)
  Sigma_samples <- rep(list(matrix(0,q,q)), niter)

  # Also list WAIC samples if return_WAIC=TRUE
  if(return_WAIC==TRUE)
    logLL_samples <- matrix(0, nrow=n, ncol=niter)
  
  # Start Gibbs sampler
  for(bj in 1:niter){
    
    if (bj %% 100 == 0)
      cat("Iteration:", bj,"/",niter, "\n")
    
    # Update nxq polya-Gamma weight matrix W
    for(j in 1:q){
      thetaj <- X%*%B[,j]+latentU[,j]
      W[,j] <- switch(response_types[j],
                      binary = pmax(BayesLogit::rpg(n,1, thetaj), nugget), # for numerical stability
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
    
    # Update latent matrix Z for Z = (Y-U)/W; W==1 if Y is continuous
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
    
    # Update latent variables matrix U and logLL if return_WAIC=TRUE
    tmp_mat <- Z-X%*%B
    for(i in 1:n){
      Omega_i <- diag(W[i,])
      mm <- Omega_i%*%tmp_mat[i,]
      JJ <- chol2inv(chol(Omega_i+Sigma))
      latentU[i,] <- mvtnorm::rmvnorm(1, JJ%*%mm, JJ)
      if(return_WAIC==TRUE){
        logLL_samples[i,bj] = -0.5*sum(2*log(diag(chol(Omega_i+Sigma)))) -0.5*t(tmp_mat[i,])%*%JJ%*%tmp_mat[i,]
      }
    }
    
    # Update zeta_i's and nu_i's
    for (i in 1:p){
      norm_term <- sum((t(B[i,]))^2) 
      v <- max(norm_term, .Machine$double.eps) # to prevent chi parameter from collapsing to 0
      zeta[i] <- GIGrvg::rgig(n=1, lambda=u-q/2, chi=v, psi=2*nu[i])
      nu[i] <- stats::rgamma(n=1, shape=a, scale=1/(tau+zeta[i]))
    }
    
    # Save samples
    B_samples[[bj]] <- B
    Sigma_samples[[bj]] <- Sigma
  }
  
  ##########################
  ## end of Gibbs sampler ##
  ##########################
  
  # Discard burn-in 
  B_samples <- utils::tail(B_samples, niter-burn) 
  Sigma_samples <- utils::tail(Sigma_samples, niter-burn) 
  
  # Calculate the WAIC if return_WAIC=TRUE
  if(return_WAIC==TRUE){ 
    logLL_samples <- logLL_samples[, (niter-burn+1):niter]
    post_mean_logLL <- rowMeans(logLL_samples)  
    post_var_logLL <- apply(logLL_samples, 1, var)
    WAIC <- -2*sum(post_mean_logLL) + 2*sum(post_var_logLL)
  }
  
  if(return_WAIC==TRUE){
    output <- list(B_samples = B_samples,
                   Sigma_samples = Sigma_samples,
                   WAIC = WAIC)
  } else {
    output <- list(B_samples = B_samples,
                   Sigma_samples = Sigma_samples)
  }
  
  # Return list
  return(output)
}


################################
## EXTRACT SUMMARY STATISTICS ##
################################

# Set threshold = 0 for one-step algorithm
# Otherwise threshold should be greater than 0
Mt_MBSP_summary = function(B_samples, Sigma_samples, threshold=0){

  p = nrow(B_samples[[1]])
  q = ncol(B_samples[[1]])
 
  # Extract the posterior median and 2.5th and 97.5th quantiles of B 
  if(p>1){
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
  } else {
    # if p=1
    B_est <- matrix(0, nrow=1, ncol=q)
    B_lower <- matrix(0, nrow=1, ncol=q)
    B_upper <- matrix(0, nrow=1, ncol=q)
    
    for(j in 1:q){
      B1j <- unlist(lapply(B_samples, "[", 1, j))
      
      # Take posterior median as point estimate for B
      B_est[1,j] <- as.numeric(stats::quantile(B1j, prob=0.5))
      # For marginal credible intervals
      B_lower[1,j] <- as.numeric(stats::quantile(B1j, prob=0.025))
      B_upper[1,j] <- as.numeric(stats::quantile(B1j, prob=0.975))
    }
  }
  
  # set A_n
  if(threshold > 0){
    B_lower <- ifelse(B_lower > -threshold*stats::sd(B_lower), threshold*stats::sd(B_lower), B_lower)
    B_upper <- ifelse(B_upper < threshold*stats::sd(B_upper), -threshold*stats::sd(B_upper), B_upper)
  }
  # for active_predictors 
  temp <- sign(B_lower)*sign(B_upper)
  active_predictors <- 1*(temp>=0)

  # Point estimates and 95% posterior credible intervals for B 
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

  # point estimates and 95% posterior credible intervals for Sigma
  Sigma_est <- matrix(Sigma_est, q, q, byrow=F)
  Sigma_lower <- matrix(Sigma_lower, q, q, byrow=F)
  Sigma_upper <- matrix(Sigma_upper, q, q, byrow=F)

  output <- list(B_est = B_est,
                 B_active = active_predictors,
                 B_lower = B_lower,
                 B_upper = B_upper,
                 Sigma_est = Sigma_est, 
                 Sigma_lower = Sigma_lower,
                 Sigma_upper = Sigma_upper)
  return(output)
}
