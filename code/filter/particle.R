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
#library(nloptr)

# ----------------------------------------------------------------------
# (Multivariate) Particle filter
# ----------------------------------------------------------------------
particle.filter <- function(y, d=ncol(data.frame(y)), var.eps=1, cov.eta=diag(d), eta.sim, u.sim, a1=0) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- nrow(u.sim)
    
    # initialize series and loglikelihood
    x.pr <- rep(0,T)
    x.up <- rep(0,T)
    loglik <- 0 
    
    # initial state estimator
    x_up <- rnorm(P, a1, 1)
    
    # estimate states (alphas) of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr <- x_up + cov.eta * eta.sim[,t] # need to simulate using cov.eta
        x.pr[t] <- mean(x_pr)
        
        # [2] (Update step)
        lik <- dnorm( y[t,]*rep(1,P) , x_pr , var.eps) 
        log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
        
        if (log.mean.lik == -Inf) {
            break
        } 
        loglik <- loglik - log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        x_up <- sir(x_pr, lik, u.sim[,t])
        x.up[t] <- mean(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, loglik=loglik))
}

# sequential importance sampling algorithm
sir <- function(x_pr_, x_wt, u) {
    
    P <- length(x_pr_)
    x_up <- rep(0,P)
    
    # sorting and weighting
    x_wt <- x_wt/sum(x_wt)
    x_sort <- cbind(seq(1,P,1),c(x_pr_)) # update for multivariate
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
    T <- length(y)
    
    # draw noise
    eta_sim <- rnorm(P*T, 0,1) # TBD mvrnorm 3-dimensional
    eta_sim <- matrix(eta_sim, nrow=P, ncol=T) 
    u_sim <- runif(P*T, 0,1)
    u_sim <- matrix(u_sim, P, T)    
    
    # set optimization parameters
    lb <- c(rep(0.1,(D-1)), -1);
    ub <- c(rep(5,  (D-1)),  1);
    theta_start <- c(rep(1,(D-1)), 0)
    obj <- function(theta){ return( -particle.filter(y, cov.eta=construct.cov(theta[1:(D-1)], theta[D]))$loglik ) } 
    
    # run box-constrained optimization
    print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- kalman.filter(y, cov.eta=construct.cov(theta_mle[1:(D-1)], theta_mle[D]))$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
