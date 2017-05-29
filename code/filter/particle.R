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
library(mvtnorm)
library(nloptr)

# ----------------------------------------------------------------------
# (Univariate) Particle filter
# ----------------------------------------------------------------------
particle.filter <- function(y, var.eps=1, cov.eta=1, eta.sim, u.sim, a1=0, P1=1) {
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
    P <- nrow(u.sim)
    
    # initialize series and loglikelihood
    x.pr <- data.frame(matrix(ncol = D, nrow = T)) 
    x.up <- data.frame(matrix(ncol = D, nrow = T)) 
    lik <- rep(0, P)
    loglik <- 0 
    
    # initial state estimator
    x_pr <- x_up <- mvrnorm(P, mu=a1, Sigma=P1)
    
    # compute Cholesky decomposition of cov.eta
    if (is.symmetric.matrix(cov.eta) && is.positive.definite(cov.eta)) {
        uL <- chol(cov.eta)
    } else {
        return(list(x.pr=NaN, x.up=NaN, loglik=-Inf))
    }

    # estimate states (alphas) of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr <- x_up + eta.sim[[t]] %*% uL
        x.pr[t,] <- colMeans(x_pr)
        
        # [2] (Update step)
        lik <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x_pr[p,], sigma=sqrt(var.eps * diag(D)))) 
        log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
 
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        x_up <- m.sir(x_pr, lik, u.sim[,t])
        x.up[t] <- mean(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, loglik=loglik))
}

# ----------------------------------------------------------------------
# Univariate sequential importance resampling algorithm
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
# Multivariate sequential importance resampling algorithm
# ----------------------------------------------------------------------
m.sir <- function(x_pr_, x_wt, u) {
    
    D <- ncol(x_pr_)
    P <- nrow(x_pr_)
    
    x_up <- matrix(0,P,D)
    
    # weighting
    x_wt <- x_wt/sum(x_wt) # normalize weights
    x_cwt <- c(0, cumsum(x_wt))
    
    j <- 1
    for (i in 1:P){
        while((x_cwt[i] < u[j]) && (u[j] <= x_cwt[i+1])){
            x_up[j,] <- x_pr_[i,]
            if (j < P) {
                j <- j+1
            }
            else break
        }
    }
    return(x_up)
}

# ----------------------------------------------------------------------
# (Univariate) maximum likelihood estimator
# ----------------------------------------------------------------------
particle.mle <- function(y, P) {
    T <- length(y)
    
    # draw noise and particles
    eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=P, ncol=T) 
    u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
    for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
    
    # set optimization parameters
    lb <- 0.1
    ub <- 5
    theta_start <- 0
    obj <- function(theta){ return( -particle.filter(y, cov.eta=theta, eta.sim=eta.sim, u.sim=u.sim)$loglik ) } 
    
    # run box-constrained optimization
    print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- particle.filter(y, cov.eta=theta_mle, eta.sim=eta.sim, u.sim=u.sim)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}

# ----------------------------------------------------------------------
# (Multivariate) maximum likelihood estimator
# ----------------------------------------------------------------------
m.particle.mle <- function(y, D, P) {
    T <- nrow(y)
    
    # draw noise and particles
    eta.sim <- list()
    for (t in 1:T) {
        eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D)) 
    }
    u.sim <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
    for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
    
    # set optimization parameters
    lb <- c(rep(0.1,D), -1)
    ub <- c(rep(5,  D),  1)
    theta_start <- c(rep(1,D), 0)
    obj <- function(theta){ return( -m.particle.filter(y, cov.eta=construct.cov(c(theta[1:D]), theta[D+1]), eta.sim=eta.sim, u.sim=u.sim)$loglik ) } 
    
    # run box-constrained optimization
    print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- m.particle.filter(y, cov.eta=construct.cov(theta_mle[1:D], theta_mle[D+1]), eta.sim=eta.sim, u.sim=u.sim)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
