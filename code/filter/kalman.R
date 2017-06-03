# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Kalman filter (for multivariate local level model)
#
# (Author) Hans-Peter Höllwirth
# (Date)   12.05.2017

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
library(nloptr)

# ----------------------------------------------------------------------
# (Multivariate) Kalman filter
# ----------------------------------------------------------------------
kalman.filter <- function(y, d=ncol(data.frame(y)), var.eps=1, cov.eta=diag(d), a1=0, P1=diag(d)) {
    y <- data.frame(y)
    T <- nrow(y)
    
    # initialize series and loglikelihood
    v <- data.frame(matrix(ncol = d, nrow = T))
    a.t <- data.frame(matrix(ncol = d, nrow = T+1)) 
    a <- data.frame(matrix(ncol = d, nrow = T)) 
    F <- K <- P <- P.t <- list()
    
    # initial state estimator
    a[1,] <- a1
    P[[1]] <- P1
    
    # estimate states (alphas) of state space model
    for (t in 1:T) {
        # [1] (Update step)
        # compare prediction to observation
        v[t,] <- y[t,] - a[t,]     
        F[[t]] <- P[[t]] + var.eps*diag(d)
        
        # compute Kalman gain (determines how much the new observation will affect the updated prediction)
        K[[t]] <- P[[t]] %*% solve(F[[t]])
        
        # compute filtered estimator to update (a posteriori) state estimate
        a.t[t,] <- a[t,] + t(K[[t]] %*% t(v[t,]))
        P.t[[t]] <- P[[t]] - K[[t]] %*% P[[t]]
        
        # [2] (Prediction step) 
        # one-step ahead (a priori) prediction
        a[t+1,] <- a.t[t,]
        P[[t+1]] <- P.t[[t]] + cov.eta
    }
    P[[T+1]] <- NULL
    
    # compute log-likelihood
    loglik <- -T*d/2 * log(2*pi)
    for (t in 1:T) {
        # TBD: add try-catch
        loglik <- loglik - 0.5 * (log(det(F[[t]])) + c(t(v[t,])) %*% solve(F[[t]]) %*% c(t(v[t,]))) 
    }
    
    return(list(a=a[1:T,], P=P, loglik=loglik))
}

# ----------------------------------------------------------------------
# Maximum likelihood estimator
# ----------------------------------------------------------------------
kalman.mle <- function(y, D, verbose=TRUE) {
    T <- length(y)
    
    # set optimization parameters
    if (D == 1) {
        lb <- 0.1
        ub <- 5
        theta_start <- 0
        obj <- function(theta){ return( -kalman.filter(y, cov.eta=theta)$loglik ) } 
        
    } else {
        lb <- c(rep(0.1,D), -1)
        ub <- c(rep(5,  D),  1)
        theta_start <- c(rep(1,D), 0)
        obj <- function(theta){ return( -kalman.filter(y, cov.eta=construct.cov(theta[1:D], theta[D+1]))$loglik ) } 
    }
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    if (D == 1) 
        loglik <- kalman.filter(y, cov.eta=theta_mle)$loglik
    else
        loglik <- kalman.filter(y, cov.eta=construct.cov(theta_mle[1:D], theta_mle[D+1]))$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
