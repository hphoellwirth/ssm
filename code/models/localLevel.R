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

# load libraries
library(MASS)


# ----------------------------------------------------------------------
# Generate univariate local level model data
# ----------------------------------------------------------------------
gen.llm.data <- function(n, var.eps=1, var.eta=1, a1=0, P1=1) {
    y <- rep(0,n)   # observations
    x <- rep(0,n+1) # states
    
    # draw intial state
    x[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    eps <- rnorm(n, mean=0, sd=sqrt(var.eps))
    eta <- rnorm(n, mean=0, sd=sqrt(var.eta))
    
    # compute observations and states
    for (t in 1:n) {
        y[t]   <- x[t] + eps[t] # observation
        x[t+1] <- x[t] + eta[t] # state
    }
    return(list(y=y, x=x[1:n]))
}


# ----------------------------------------------------------------------
# Generate multivariate local level model data
# ----------------------------------------------------------------------
gen.multi.llm.data <- function(n, d=2, var.eps=1, cov.eta=diag(d), a1=0, P1=diag(d)) {
    y <- data.frame(matrix(ncol = d, nrow = n)) # observations
    x <- data.frame(matrix(ncol = d, nrow = n)) # states
    
    # draw intial state
    x[1,] <- mvrnorm(1, mu=rep(a1,d), Sigma=P1)
    
    # draw disturbances
    eps <- mvrnorm(n, mu=rep(0,d), Sigma=var.eps*diag(d))
    eta <- mvrnorm(n, mu=rep(0,d), Sigma=cov.eta) 
    
    # compute observations and states
    for (t in 1:n) {
        y[t,]   <- x[t,] + eps[t,] # observation
        x[t+1,] <- x[t,] + eta[t,] # state
    }
    return(list(y=y, x=x[1:n,]))
}



