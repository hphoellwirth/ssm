# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Kalman filter
#
# (Author) Hans-Peter Höllwirth
# (Date)   13.04.2017


# ----------------------------------------------------------------------
# (Multivariate) Kalman filter
# ----------------------------------------------------------------------
kalman.filter <- function(y, d=ncol(data.frame(y)), var.eps=1, cov.eta=diag(d), a1=0, P1=diag(d)) {
    y <- data.frame(y)
    n <- nrow(y)
    
    # initialize series and loglikelihood
    v <- data.frame(matrix(ncol = d, nrow = n))
    a.t <- data.frame(matrix(ncol = d, nrow = n+1)) 
    a <- data.frame(matrix(ncol = d, nrow = n)) 
    F <- K <- P <- P.t <- list()
    
    l <- rep(0,n)
    loglik <- 0
    
    # initial state estimator
    a[1,] <- a1
    P[[1]] <- P1
    
    # estimate states (alphas) of state space model
    for (t in 1:n) {
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
        
        # [3] (Log-likelihood computation)
        loglik <- loglik - 0.5 * (d*log(2*pi) + log(det(F[[t]])) + c(t(v[t,])) %*% solve(F[[t]]) %*% c(t(v[t,])))
        l[t] <- loglik
    }
    P[[n+1]] <- NULL
    return(list(a=a[1:n,], P=P, l=l))
}

