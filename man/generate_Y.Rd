% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generation.R
\name{generate_Y}
\alias{generate_Y}
\title{Generates n-by-q response matrix according to the settings in Section 5.1 of the paper. The responses  are in the q columns and can be of mixed-type (continuous, binary, or count).}
\usage{
  generate_Y(X, B0, Sigma0, response_types, r=50)
}
\arguments{
  \item{X}{n-by-p design matrix of covariates}
  
  \item{B0}{p-by-q ground truth regression coefficients matrix}
  
  \item{Sigma0}{q-by-q ground truth covariance matrix}
  
  \item{response_types}{vector of q response types. The entries of the vector must be one of the following: "continuous", "binary", or "count".}
  
  \item{r}{number of successful trials. Only for count responses. Default is r=50 for equidispersed data. To obtain overdispersed data, decrease r.}
  
}

\value{An n-by-q response matrix Y. Each column of Y corresponds to the responses of the same type (i.e. continuous, binary, or count).
}

\description{
This function generates an n-by-q response matrix of mixed-type responses. Continuous responses are generated according to a Gaussian distribution, binary responses are generated according to a Bernoulli distribution, and count responses are generated according to a negative binomial distribution. The simulation settings are those given in Section 5.1 of the manuscript. 
}

\examples{
n <- 100  # sample size 
p <- 20   # number of covariates
s <- 5    # number of significant covariates
q <- 3    # number of multiple responses 
response_types <- c('binary','continuous','count')

# Set seed to reproduce results later
set.seed(1)

# Generate p-by-q regression coefficients matrix B0, where s=# of nonzero rows.
B0 <- generate_B0(p, s, response_types)
B0

# Generate true q-by-q covariance matrix Sigma0
Sigma0 <- generate_Sigma0(q)
# Generate nxp design matrix X
X <- generate_X(n, p)

# Finally generate n-by-q response matrix Y
Y <- generate_Y(X, B0, Sigma0, response_types)
}