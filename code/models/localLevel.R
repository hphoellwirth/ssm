# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Local level models
#
# (Author) Hans-Peter Höllwirth
# (Date)   28.04.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# house cleaning
rm(list = ls())
par(mfrow=c(1,1))

# load libraries

# if interactive, during the development, set to TRUE
interactive <- TRUE
if (interactive) {
    setwd("/Users/Hans-Peter/Documents/Masters/14D000/code/models")
} 


# ----------------------------------------------------------------------
# Generate univariate local level model data
# ----------------------------------------------------------------------
gen.llm.data <- function(n, var.eps=1, var.eta=1, a1=0, P1=1) {
    y <- rep(0,n)
    alpha <- rep(0,n+1)
    
    # draw intial state
    alpha[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    eps <- rnorm(n, mean=0, sd=sqrt(var.eps))
    eta <- rnorm(n, mean=0, sd=sqrt(var.eta))
    
    # compute observations and states
    for (t in 1:n) {
        y[t]       <- alpha[t] + eps[t] # observation
        alpha[t+1] <- alpha[t] + eta[t] # state
    }
    return(list(y=y, alpha=alpha[1:n]))
}

# ----------------------------------------------------------------------
# Generate trivariate local level model data
# ----------------------------------------------------------------------
gen.tllm.data <- function(n, var.eps=1, var.eta=1, a1=0, P1=1) {
    y <- data.frame(matrix(ncol = 3, nrow = n))
    alpha <- data.frame(matrix(ncol = 3, nrow = n))
    
    # draw intial state
    alpha[1,] <- rnorm(3, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    eps <- rnorm(3*n, mean=0, sd=sqrt(var.eps))
    eta <- rnorm(3*n, mean=0, sd=sqrt(var.eta)) # TBD: covariance matrix
    
    # compute observations and states
    for (t in 1:n) {
        y[t,]       <- alpha[t,] + eps[t:(t+2)] # observation
        alpha[t+1,] <- alpha[t,] + eta[t:(t+2)] # state
    }
    return(list(y=y, alpha=alpha[1:n]))
}



