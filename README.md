# MtMBSP

This is an R package to implement the Mixed-typed Multivariate Bayesian Model with Shrinkage Priors (Mt-MBSP) introduced by Shao-Hsuan Wang, Ray Bai, and Hsin-Hsiung Huang in their paper "Mixed-type Multivariate Bayesian Sparse Variable Selection with Shrinkage Priors" (pre-print, 2023+).

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
        u=0.5, a=0.5, tau=.001, d1 = dim(Y)[2], d2=10, c1=10, c2=1,
        algorithm = "1step",
        step1_iter = 200, bound_error=0, 
        max_iter = 2000, burnin=1000, details=TRUE)
```

`X` is the n-by-p design matrix with n samples of p covariates.

`Y` is the n-by-q response matrix with n samples of q mixed-type responses. 

`response_types` is a q-dimensional string vector which indicates the different response types. This vector must only contain entries of `"continuous"`, `"binary"`, or `"count"`. 
Continuous responses are modeled with the Gaussian distribution, binary responses are modeled with the Bernoulli distribution, and count responses are modeled with the negative binomial distribution.
It is extremely important that this vector give the correct response types corresponding to the q columns of `Y`; otherwise, the code may fail to run.

For example, if `Y` contains a count variable in its first column, a binary variable in its second column, a continuous variable in the third column, another continuous variable in the fourth column, and another binary variable 
in its fifth column, then we should pass `response_types = c("count","binary","continuous","continuous","binary")` to the `Mt_MBSP` function.

`u` and `a` are the hyperparameters in the TPBN prior for the local scale parameters. `u=0.5` and `a=0.5` correspond to the popular horseshoe prior of Carvalho et al. (2010).  Meanwhile, `tau` is the global hyperparameter in the TPBN prior.

`d1` and `d2` are the degrees of freedom and the scale parameter in the inverse-Wishart prior on Sigma.

`c1` and `c2` are the shape and rate parameters for the gamma prior on the unknown dispersion parameters `r` for the negative binomial responses. These arguments are ignored for non-count responses.

`algorithm` indicates whether to use the one-step approach (`"1step"`) or two-step approach (`"2step"`) for variable selection. The default is `"1step"`. However, for large p, it may be advantageous
to use the two-step approach. If it is desired to use the two-step algorithm, then this argument should be `algorithm="2step"`.

`step1_iter` indicates the number of iterations to run Step 1 of the two-step approach. This argument is only used if `algorithm="2step"`, and it must be less than `max_iter`.

`bound_error` is the threshold used to screen variables in Step 1 of the two-step approach. If `bound_error=0` is specified, then there is no variable screening step done.
If `algorithm="2step"`, then the `Mt_MBSP` function will use a default of `bound_error=0.02` in Step 1. However, the user can specify a different nonnegative `bound_error` argument. Note, however, that `bound_error` is only used if `algorithm="2step"`.

`max_iter` is the total number of Gibbs sampling iterations to run.

`burnin` is the number of burn-in samples, so that the final total number of MCMC samples to approximate the posterior is `max_iter-burnin`.

`details` indicates whether to return all of the saved posterior samples and sample quantiles. By default, `details=TRUE`.

## 3. Additional Functions for Simulating Synthetic Data

In addition, there are several functions `generate_B0`, `generate_Sigma0`, `generate_X`, and `generate_Y`, which are used to generate synthetic datasets for the simulation studies in Section 5.1 of the manuscript.

We demonstrate the usage of these functions, in addition to the main `Mt_MBSP` function in the simulated examples below.

## 4. Example of the One-step Approach

To demonstrate the `Mt_MBSP` function, we first illustrate it on a small synthetic dataset with n=100, p=20, s=5 (where s is the number of significant covariates), and q=3 responses which are binary, continuous, and count, i.e. `response_types = c('binary','continuous','count')`.

```
n <- 100  # sample size 
p <- 20   # number of covariates
s <- 5    # number of significant covariates
q <- 3    # number of multiple responses 
response_types <- c('binary','continuous','count')

# Set seed to reproduce results later
set.seed(123)

# Generate p-by-q regression coefficients matrix B0, where s=# of nonzero rows.
B0 <- generate_B0(p, s, response_types)
# Generate true q-by-q covariance matrix Sigma0
Sigma0 <- generate_Sigma0(q)
# Generate n-by-p design matrix X
X <- generate_X(n, p)
# Generate n-by-q response matrix Y
Y <- generate_Y(X, B0, Sigma0, response_types)
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
# rMSE = 0.172689

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
classifications <- rep(0, p)
selected_variables <- which(rowSums(output$B_active)!=0)
classifications[selected_variables] <- 1

# Ground truth
truth <- rep(0, p)
true_nonzero_variables <- which(rowSums(B0)!=0)
truth[true_nonzero_variables] <- 1

# Compare selected variables to the ground truth significant variables 
selected_variables 
# 3 10 14 15 19
true_nonzero_variables
# 3 10 14 15 19
```

## 5. Example of the Two-step Approach for Large p

When p is large, it may be advantageous to use the two-step approach. The first step of our two-stage algorithm screens out a large number of variables using only a tiny number of Gibbs sampling steps and a thresholding rule determined by the argument `bound_error`.
In the second step, the remaining regression coefficients are estimated with the Gibbs sampler.

The computational gains from using two-step algorithm rather than one-step are much more evident for larger p. However, the two-step estimator often has superior estimation and variable selection over the one-step estimator when p is large, even if the computational gain is modest.
Below, we demonstrate how to use the `Mt_MBSP` function with the two-step algorithm on a simulated dataset with n=150, p=1000, s=10 (where s is the number of significant covariates), and q=4 responses with two continuous and two binary responses, i.e. `response_types = c('continuous','binary','continuous','binary')`.

```
n <- 150  # sample size 
p <- 1000 # number of covariates
s <- 10   # number of significant covariates
q <- 4    # number of multiple responses 

# Responses consist of 2 continuous and 2 binary outcomes  
response_types <- c('continuous','binary','continuous','binary')

# Set seed to reproduce results later
set.seed(123)

# Generate p-by-q regression coefficients matrix B0, where s=# of nonzero rows.
B0 <- generate_B0(p, s, response_types)
# Generate true q-by-q covariance matrix Sigma0
Sigma0 <- generate_Sigma0(q)
# Generate n-by-p design matrix X
X <- generate_X(n, p)
# Generate n-by-q response matrix Y
Y <- generate_Y(X, B0, Sigma0, response_types)
```

We next fit the two-step algorithm. In this case, we must specify `algorithm='2step'` in the `Mt_MBSP` function. We run Step 1 of the two-step approach with 200 iterations (i.e. `step1_iter=200`) and a threshold `bound_error=0.02`. 
The total number of iterations is `max_iter=2000`, which indicates that we run Step 2 for 1800 iterations.

```
# Fit two-step model. In this case, you should specify the following arguments:
#    algorithm='2step' 
#    step1_iter (which should be less than burnin), 
#    bound_error (threshold gamma)

response_types <- c('continuous','binary','continuous','binary')

output <- Mt_MBSP(X, Y, response_types,
                  algorithm='2step', step1_iter=200, bound_error=0.02,
                  max_iter=2000, burnin=1000)

# VERY IMPORTANT: Need to make sure that the response_types argument in Mt_MBSP is correctly specified, 
# otherwise the code may not run correctly!
```

The main results of interest might be: `output$B_est` (posterior median estimate for regression coefficients matrix B), `output$B_active` (a binary matrix with `1` indicating that the variable is selected, and `0` indicating that it is not selected), `output$B_lower` (the 0.025 quantiles for the entries of B), and `output$B_upper` (the 0.975 quantiles for the entries of B).

We can also obtain the following performance metrics.

```
# root mean squared error (rMSE) for two-step estimator
rMSE <- sqrt(sum((output$B_est-B0)^2)/(p*q))
rMSE
# rMSE = 0.06893186

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
# CP = 0.9995

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
# 118 179 195 229 299 415 463 526 818 938
true_nonzero_variables
# 118 179 195 229 299 415 463 526 818 938
```

