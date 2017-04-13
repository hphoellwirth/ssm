# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Kalman filter
#
# (Author) Hans-Peter Höllwirth
# (Date)   13.04.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# house cleaning
#rm(list = ls())
par(mfrow=c(1,1))

# load libraries

# if interactive, during the development, set to TRUE
interactive <- TRUE
if (interactive) {
    setwd("/Users/Hans-Peter/Documents/Masters/14D000/code/filter")
} 


# ----------------------------------------------------------------------
# Kalman filter
# ----------------------------------------------------------------------
kalman.filter <- function(y, var.eps=1, var.eta=1, a1=0, P1=1) {
    n <- length(y)
    
    # initialize series
    v <- F <- K <- rep(0,n)
    a.t <- P.t <- rep(0,n)
    a <- P <- rep(0,n+1)
    
    # initial state estimator
    a[1] <- a1
    P[1] <- P1
    
    # estimate states (alphas) of state space model
    for (t in 1:n) {
        # one-step ahead prediction error of y
        v[t] <- y[t] - a[t]
        F[t] <- P[t] + var.eps
        
        # Kalman gain
        K[t] <- P[t] / F[t]
        
        # filtered estimator
        a.t[t] <- a[t] + K[t]*v[t]
        P.t[t] <- P[t] * (1 - K[t])
        
        # one-step ahead prediction
        a[t+1] <- a.t[t]
        P[t+1] <- P.t[t] + var.eta
        
    }
    return(list(a=a[1:n], P=P[1:n]))
}

