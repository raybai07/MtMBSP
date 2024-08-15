# MtMBSP

This is an R package to implement the Mixed-typed Multivariate Bayesian Model with Shrinkage Priors (Mt-MBSP) introduced by Shao-Hsuan Wang, Ray Bai, and Hsin-Hsiung Huang in their paper "Two-Step Mixed-Type Multivariate Bayesian Sparse Variable Selection with Shrinkage Priors" (preprint: https://arxiv.org/abs/2201.12839).

Our model implements Bayesian variable selection and estimation for mixed-type (i.e. mixture of continuous and discrete) multivariate regression using shrinkage priors. Specifically, the three parameter beta normal (TPBN) family of priors is used, which contains the horseshoe prior (Carvalho et al., 2010) as a special case. Our method enables joint analysis of mixed outcomes and facilitates variable selection when the number of covariates can be much larger than sample size. 
For large p, we introduce a two-step variable selection approach that can also be implemented with this R package.

## 1. Installing the Package

In order to get this code to run, you must have R (>= 3.6.0) installed. It would be ideal to have the most recent version of R. You can use the following R code to install and load the MtMBSP package. 

```
# Install MtMBSP package
library(devtools)
install_github(repo = "raybai07/MtMBSP")

# Load MtMBSP package
library(MtMBSP)
``` 
## 2. Main Function

The main function to implement our method is `Mt_MBSP` which we will describe in detail and demonstrate how to use below. 

```
Mt_MBSP(X, Y, response_types,
        u=0.5, a=0.5, tau=1/(dim(X)[2]*sqrt(dim(X)[1]*log(dim(X)[1]))),
        d1 = dim(Y)[2], d2=10, c1=10, c2=1,
        algorithm = c("1step", "2step"),
        niter = 1100, burn=100,
        step2_niter = 1100, step2_burn=100, 
        threshold = seq(from=0.02, to=0.40, by=0.02),
        parallelize = TRUE, ncores = 10)
```

`X` is the n-by-p design matrix with n samples of p covariates.

`Y` is the n-by-q response matrix with n samples of q mixed-type responses. 

`response_types` is a q-dimensional string vector which indicates the different response types. This vector must only contain entries of `"continuous"`, `"binary"`, or `"count"`. 
Continuous responses are modeled with the Gaussian distribution, binary responses are modeled with the Bernoulli distribution, and count responses are modeled with the negative binomial distribution.
It is extremely important that this vector give the correct response types corresponding to the q columns of `Y`; otherwise, the code may fail to run.

For example, if `Y` contains a count variable in its first column, a binary variable in its second column, a continuous variable in the third column, another continuous variable in the fourth column, and another binary variable 
in its fifth column, then we should pass `response_types = c("count","binary","continuous","continuous","binary")` to the `Mt_MBSP` function.

`u` and `a` are the hyperparameters in the TPBN prior for the local scale parameters. `u=0.5` and `a=0.5` correspond to the popular horseshoe prior of Carvalho et al. (2010).  Meanwhile, `tau` is the global hyperparameter in the TPBN prior. The default for `tau` is 1/(p * sqrt(n * log(n))). If this quantity becomes lower than `1e-5`, then `tau` is set as `1e-5` for numerical stability reasons. If p is much than n and all of the signals in the data are very weak, then our method may select a null model. In this case, it is recommended that the user increase the value of `tau`.

`d1` and `d2` are the degrees of freedom and the scale parameter in the inverse-Wishart prior on Sigma.

`c1` and `c2` are the shape and rate parameters for the gamma prior on the unknown dispersion parameters `r` for the negative binomial responses. These arguments are ignored for non-count responses.

`algorithm` indicates whether to use the one-step approach (`"1step"`) or two-step approach (`"2step"`) for variable selection. The default is `"1step"`. However, for large p, it may be advantageous
to use the two-step approach. If it is desired to use the two-step algorithm, then this argument should be `algorithm="2step"`.

`niter` indicates the number of MCMC iterations to run in the one-step algorithm OR in Step 1 of the two-step algorithm. `burn` is the number of MCMC samples to discard as burnin. 

`step2_niter` indicates the number of MCMC iterations to run in Step 2 of the two-step algorithm, and `burn` is the number of MCMC to discard as burnin in Step 2. These arguments are only used if `algorithm="2step"`.

`threshold` is a grid of thresholds gamma to search over in Step 1 of the two-step algorithm. The threshold which minimizes the WAIC is used as the final model in Step 2. The default is an equispaced grid 0.02, 0.04, ..., 0.40. This argument is only used if `algorithm="2step"`.

`parallelize` is a Boolean variable for whether to parallelize Step 2 of the two-step algorithm over the grid of values in `threshold`. `ncores` is the number of cores to use for parallelization. These arguments are only used if `algorithm="2step"`. 


## 3. Additional Functions for Simulating Synthetic Data

In addition, there are several functions `generate_B0`, `generate_X`, and `generate_Y`, which are used to generate synthetic datasets using the settings in Section 5.1 of the paper.

We demonstrate the usage of these functions, in addition to the main `Mt_MBSP` function in the simulated examples below.

## 4. Example of the One-step Approach

To demonstrate the `Mt_MBSP` function, we first illustrate it on a small synthetic dataset with n=100, p=20, s=5 (where s is the number of significant covariates), and q=3 responses which are binary, continuous, and count, i.e. `response_types = c('binary','continuous','count')`.

```
n <- 100  # sample size 
p <- 500  # number of covariates
s <- 5    # number of significant covariates
q <- 3    # number of multiple responses 
response_types <- c('binary','continuous','count')

# Set seed to reproduce results later
set.seed(1234)

# Generate p-by-q regression coefficients matrix B0, where s=# of nonzero rows.
B0 <- generate_B0(p, s, response_types)
# Generate n-by-p design matrix X
X <- generate_X(n, p)
# Generate n-by-q response matrix Y
Y <- generate_Y(X, B0, response_types)
```

We then fit the Mt-MBSP model as follows using the main `Mt_MBSP` function.

```
# Fit one-step model with X, Y, and response_types, and use default arguments for all other arguments.

response_types <- c('binary','continuous','count')

output <- Mt_MBSP(X, Y, response_types)

# VERY IMPORTANT: Need to make sure that the response_types argument in the Mt_MBSP function  
# is correctly specified, otherwise the code may not run correctly!
```

The main results of interest might be: `output$B_est` (posterior median estimate for regression coefficients matrix B), `output$B_active` (a binary matrix with `1` indicating that the variable is selected, and `0` indicating that it is not selected), `output$B_lower` (the 0.025 quantiles for the entries of B), and `output$B_upper` (the 0.975 quantiles for the entries of B).

We can also obtain the following performance metrics.

```
# root mean squared error (rMSE) for one-step estimator
rMSE <- sqrt(sum((output$B_est-B0)^2)/(p*q))
rMSE
# rMSE = 0.03186077

# Coverage probability (CP) for one-step estimator
coverage_mat <- matrix(0, nrow=p, ncol=q)
for(j in 1:p){
    for(k in 1:q){
        if(B0[j,k]>=output$B_lower[j,k] & B0[j,k]<=output$B_upper[j,k])
           coverage_mat[j,k] <- 1
    }
}
CP <- sum(coverage_mat)/(p*q)
CP  
# CP = 1

# Variable selection performance for one-step estimator
classifications_mat <- output$B_active
 
TP <- as.numeric(length(which(classifications_mat==1 & B0!=0)))
TN <- as.numeric(length(which(classifications_mat==0 & B0==0)))
FP <- as.numeric(length(which(classifications_mat==1 & B0==0)))
FN <- as.numeric(length(which(classifications_mat==0 & B0!=0)))
 
sens <- TP/(TP+FN)
spec <- TN/(TN+FP)
MCC <- (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))

sens
# sens = 0.8181818
spec
# spec = 1
MCC
# MCC = 0.9039272
```

## 5. Example of the Two-step Approach for Large p

When p is large, it may be advantageous to use the two-step approach. The first step of our two-stage algorithm screens out a large number of variables using Gibbs sampling and a thresholding rule. An optimal threshold is chosen from a grid of candidate values in `threshold'. In the second step, the remaining regression coefficients are estimated with the Gibbs sampler. 

Note that if p is much larger than n and the signals are all very weak, then the two-stage estimator may return a null model. In this case, it is recommended that the user increase the maximum candidate threshold in `threshold'.

Below, we demonstrate how to use the `Mt_MBSP` function with the two-step algorithm on a simulated dataset with n=150, p=1000, s=10 (where s is the number of significant covariates), and q=4 responses with two continuous and two binary responses, i.e. `response_types = c('continuous','binary','continuous','binary')`.

```
n <- 150  # sample size 
p <- 1000 # number of covariates
s <- 10   # number of significant covariates
q <- 4    # number of multiple responses 

# Responses consist of 2 continuous and 2 binary outcomes  
response_types <- c('continuous','binary','continuous','binary')

# Set seed to reproduce results later
set.seed(1234)

# Generate p-by-q regression coefficients matrix B0, where s=# of nonzero rows.
B0 <- generate_B0(p, s, response_types)
# Generate true q-by-q covariance matrix Sigma0
Sigma0 <- generate_Sigma0(q)
# Generate n-by-p design matrix X
X <- generate_X(n, p)
# Generate n-by-q response matrix Y
Y <- generate_Y(X, B0, Sigma0, response_types)
```

We next fit the two-step algorithm. In this case, we must specify `algorithm='2step'` in the `Mt_MBSP` function. The argument `threshold` is a grid of candidate thresholds gamma which are tuned in Step 1 of the two-step algorithm. The value in `threshold` which minimizes the Watanabe–Akaike information criterion (WAIC) is used to select the final model in the two-step algorithm. If the two-step algorithm is used, then it may also be advantageous to parallelize Step 2 of the algorithm (`parallelize=TRUE`). In this case, the Step 2 model for each candidate value in `threshold` is computed in parallel. 

```
# Fit two-step model. In this case, you should specify algorithm='2step' 

response_types <- c('continuous','binary','continuous','binary')

ncores <- floor(parallel::detectCores()*0.75)

output <- Mt_MBSP(X, Y, response_types,
                  algorithm='2step',
                  threshold=seq(from=0.02, to=0.40, by=0.02),
		  parallelize = TRUE, ncores = ncores)

# VERY IMPORTANT: Need to make sure that the response_types argument in Mt_MBSP is correctly specified,
# otherwise the code may not run correctly!
```

The main results of interest might be: `output$B_est` (posterior median estimate for regression coefficients matrix B), `output$B_active` (a binary matrix with `1` indicating that the variable is selected, and `0` indicating that it is not selected), `output$B_lower` (the 0.025 quantiles for the entries of B), and `output$B_upper` (the 0.975 quantiles for the entries of B).

We can also obtain the following performance metrics.

```
# root mean squared error (rMSE) for two-step estimator
rMSE <- sqrt(sum((output$B_est-B0)^2)/(p*q))
rMSE
# rMSE = 0.02479676

# Coverage probability (CP) for two-step estimator
coverage_mat <- matrix(0, nrow=p, ncol=q)

for(j in 1:p){
  for(k in 1:q){
    if(B0[j,k]>=output$B_lower[j,k] & B0[j,k]<=output$B_upper[j,k])
      coverage_mat[j,k] <- 1
  }
}
CP <- sum(coverage_mat)/(p*q)
CP 
# CP = 0.99975

# Variable selection performance for two-step estimator
classifications <- rep(0, p)
selected_variables <- which(rowSums(output$B_active)!=0)
classifications[selected_variables] <- 1

# Ground truth
truth <- rep(0, p)
true_nonzero_variables <- which(rowSums(B0)!=0)
truth[true_nonzero_variables] <- 1

# Compare selected variables to the ground truth significant variables
selected_variables
# 101 284 400 623 645 848 900 905 918 934
true_nonzero_variables
# 101 284 400 623 645 848 900 905 918 934
```
