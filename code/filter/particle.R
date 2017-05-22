# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Particle filter (for multivariate local level model)
#
# (Author) Hans-Peter Höllwirth
# (Date)   18.05.2017

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
library(matrixcalc)
library(nloptr)

# ----------------------------------------------------------------------
# (Univariate) Particle filter
# ----------------------------------------------------------------------
particle.filter <- function(y, d=ncol(data.frame(y)), var.eps=1, cov.eta=diag(d), eta.sim, u.sim, a1=0, P1=1) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- nrow(u.sim)
    
    # initialize series and loglikelihood
    x.pr <- rep(0,T)
    x.up <- rep(0,T)
    loglik <- 0 
    
    # initial state estimator
    x_up <- rnorm(P, mean=a1, sd=P1)
    
    # estimate states (alphas) of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr <- x_up + sqrt(cov.eta) * eta.sim[,t] 
        x.pr[t] <- mean(x_pr)
        
        # [2] (Update step)
        # compute importance weights
        lik <- dnorm( y[t,]*rep(1,P), mean=x_pr, sd=sqrt(var.eps)) 
        log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
        
        # update likelihood
        if (log.mean.lik == -Inf) break
        loglik <- loglik + log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        x_up <- sir(x_pr, lik, u.sim[,t])
        x.up[t] <- mean(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, loglik=loglik))
}

# ----------------------------------------------------------------------
# (Multivariate) Particle filter
# ----------------------------------------------------------------------
m.particle.filter <- function(y, D=ncol(data.frame(y)), var.eps=1, cov.eta=diag(D), eta.sim, u.sim, a1=rep(0,D), P1=diag(D)) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- nrow(u.sim[[1]])
    
    # initialize series and loglikelihood
    x.pr <- data.frame(matrix(ncol = D, nrow = T)) 
    x.up <- data.frame(matrix(ncol = D, nrow = T)) 
    lik <- matrix(ncol = D, nrow = P)
    loglik <- 0 
    
    # initial state estimator
    x_up <- mvrnorm(P, mu=a1, Sigma=P1)
    
    # estimate states (alphas) of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        if (is.positive.definite(cov.eta)) {
            x_pr <- x_up +  eta.sim[[t]] %*% chol(cov.eta) # TBD: handle error where chol cannot be computed
            x.pr[t,] <- colMeans(x_pr)
            
            # [2] (Update step)
            # compute importance weights
            for (d in 1:D)
                lik[,d] <- dnorm( y[t,d]*rep(1,P) , mean=x_pr[,d] , sd=sqrt(var.eps)) 
            log.mean.lik <- tryCatch(sum(log(colMeans(lik))), error=function(e)(Inf))
        } else {
            log.mean.lik <- Inf
        }
        
        # update likelihood
        if (log.mean.lik == Inf) {
            loglik <- Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        for (d in 1:D)
            x_up[,d] <- sir(x_pr[,d], lik[,d], u.sim[[t]][,d])
        x.up[t,] <- colMeans(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, loglik=loglik))
}

# ----------------------------------------------------------------------
# Sequential importance resampling algorithm
# ----------------------------------------------------------------------
sir <- function(x_pr_, x_wt, u) {
    
    P <- length(x_pr_)
    x_up <- rep(0,P)
    
    # sorting and weighting
    x_wt <- x_wt/sum(x_wt) # normalize weights
    x_sort <- cbind(seq(1,P,1),c(x_pr_)) 
    x_pr_idx <- x_sort[order(x_sort[,2]),]
    x_pr_ <- x_pr_idx[,2]
    x_idx <- x_pr_idx[,1]
    x_wt <- x_wt[x_idx]
    x_cwt <- c(0, cumsum(x_wt))
    
    j <- 1
    for (i in 1:P){
        while((x_cwt[i] < u[j]) && (u[j] <= x_cwt[i+1])){
            x_up[j] <- x_pr_[i]
            if (j < P){
                j <- j+1
            }
            else break
        }
    }
    return(x_up)
}

# ----------------------------------------------------------------------
# Maximum likelihood estimator
# ----------------------------------------------------------------------
particle.mle <- function(y, D, P) {
    T <- nrow(y)
    
    # draw noise and particles
    eta.sim <- u.sim <- list()
    for (t in 1:T) {
        eta.sim[[t]] <- mvrnorm(P, mu=rep(0,(D-1)), Sigma=diag(D-1)) 
        u.sim[[t]]   <- matrix(runif(P*(D-1), min=0, max=1), nrow=P, ncol=(D-1)) 
        for (d in 1:(D-1)) {u.sim[[t]][,d] <- sort( u.sim[[t]][,d] )}
    }
    
    #eta_sim <- rnorm(P*T, 0,1) # TBD mvrnorm 3-dimensional
    #eta_sim <- matrix(eta_sim, nrow=P, ncol=T) 
    #u_sim <- runif(P*T, 0,1)
    #u_sim <- matrix(u_sim, P, T)    
    
    # set optimization parameters
    lb <- c(rep(0.1,(D-1)), -1);
    ub <- c(rep(5,  (D-1)),  1);
    theta_start <- c(rep(1,(D-1)), 0)
    obj <- function(theta){ return( -m.particle.filter(y, cov.eta=construct.cov(theta[1:(D-1)], theta[D]), eta.sim=eta.sim, u.sim=u.sim)$loglik ) } 
    
    # run box-constrained optimization
    print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- m.particle.filter(y, cov.eta=construct.cov(theta_mle[1:(D-1)], theta_mle[D]), eta.sim=eta.sim, u.sim=u.sim)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
