# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Local level models
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
library(MASS)


# ----------------------------------------------------------------------
# Generate univariate local level model realization
# ----------------------------------------------------------------------
gen.llm.data <- function(T, var.eps=1, var.eta=1, a1=0, P1=1) {
    y <- rep(0,T)   # observations
    x <- rep(0,T+1) # states
    
    # draw intial state
    x[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    eps <- rnorm(T, mean=0, sd=sqrt(var.eps))
    eta <- rnorm(T, mean=0, sd=sqrt(var.eta))
    
    # compute observations and states
    for (t in 1:T) {
        y[t]   <- x[t] + eps[t] # observation
        x[t+1] <- x[t] + eta[t] # state
    }
    return(list(y=y, x=x[1:T]))
}


# ----------------------------------------------------------------------
# Generate multivariate local level model realization
# ----------------------------------------------------------------------
gen.multi.llm.data <- function(T, d=2, var.eps=1, cov.eta=diag(d), a1=0, P1=diag(d)) {
    y <- data.frame(matrix(ncol = d, nrow = T)) # observations
    x <- data.frame(matrix(ncol = d, nrow = T)) # states
    
    # draw intial state
    x[1,] <- mvrnorm(1, mu=rep(a1,d), Sigma=P1)
    
    # draw disturbances
    eps <- mvrnorm(T, mu=rep(0,d), Sigma=var.eps*diag(d))
    eta <- mvrnorm(T, mu=rep(0,d), Sigma=cov.eta) 
    
    # compute observations and states
    for (t in 1:T) {
        y[t,]   <- x[t,] + eps[t,] # observation
        x[t+1,] <- x[t,] + eta[t,] # state
    }
    return(list(y=y, x=x[1:T,]))
}

# ----------------------------------------------------------------------
# Construct covariance matrix from variances and correlation 
# ----------------------------------------------------------------------
construct.cov <- function(var, rho) {
    d <- length(var)
    cov <- matrix(nrow=d, ncol=d)
    
    for (i in 1:d) {
        for (j in 1:d) {
            if (i==j)
                cov[i,j] <- var[i]
            else
                cov[i,j] <- rho * (sqrt(var[i]) * sqrt(var[j]))
        }
    }
    return(cov)
}


