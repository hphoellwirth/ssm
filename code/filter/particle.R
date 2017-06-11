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

# load samplers
source("sampler/sir.R")

# ----------------------------------------------------------------------
# (Univariate) Particle filter
# ----------------------------------------------------------------------
#y=llm.data$y; cov.eta=var.eta; x_up.init=rep(0,P)
#var.eps=1; eta.sim=NA; u.sim=NA; use.csir=TRUE

particle.filter <- function(y, var.eps=1, cov.eta=1, eta.sim=NA, u.sim=NA, P=ncol(u.sim), a1=0, P1=1, x_up.init=rnorm(P, mean=a1, sd=P1), use.csir=TRUE) {
    y <- data.frame(y)
    T <- nrow(y)
    
    # draw standard normal transition noise 
    if (is.na(eta.sim)[1][1])
        eta.sim <- matrix(rnorm(T*P, mean=0, sd=1), nrow=T, ncol=P) 
    # draw particles from standard uniform 
    if (is.na(u.sim)[1][1]) {
        u.sim   <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
        for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
    }
    
    # initialize series and loglikelihood
    x.pr <- x.up <- rep(0,T)
    x.pr.particles <- x.up.particles <- matrix(nrow = T, ncol = P)
    loglik <- 0 
    
    # initial state estimator
    x_up <- x_up.init
    
    # estimate states of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr <- x_up + sqrt(cov.eta) * eta.sim[t,] 
        x.pr.particles[t,] <- x_pr
        x.pr[t] <- mean(x_pr)
        
        # [2] (Update step)
        # compute importance weights
        lik <- dnorm( y[t,]*rep(1,P), mean=x_pr, sd=sqrt(var.eps)) 
        log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
        
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        if (use.csir) {
            x_up <- csir(x_pr, lik, u.sim[t,]) # cannot be used in combination with auxiliary filter
        } else {
            x_up <- sir(x_pr, lik, u.sim[t,])
        }

        x.up.particles[t,] <- x_up
        x.up[t] <- mean(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, x.pr.particles=x.pr.particles, x.up.particles=x.up.particles, loglik=loglik))
}

# ----------------------------------------------------------------------
# (Multivariate) Particle filter
# ----------------------------------------------------------------------
m.particle.filter <- function(y, D=ncol(data.frame(y)), var.eps=1, cov.eta=diag(D), eta.sim, u.sim, a1=rep(0,D), P1=diag(D), x_up.init=mvrnorm(P, mu=a1, Sigma=P1)) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- ncol(u.sim)
    
    # initialize series and loglikelihood
    x.pr <- x.up <- data.frame(matrix(ncol = D, nrow = T)) 
    x.pr.particles <- x.up.particles <- list()
    loglik <- 0 

    # initial state estimator
    x_up <- x_up.init
    
    # compute Cholesky decomposition of cov.eta
    if (is.symmetric.matrix(cov.eta) && is.positive.definite(cov.eta)) {
        uL <- chol(cov.eta)
    } else {
        return(list(x.pr=NaN, x.up=NaN, loglik=-Inf))
    }

    # estimate states of state space model
    for (t in 1:T) {
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr <- x_up + eta.sim[[t]] %*% uL
        x.pr.particles[[t]] <- x_pr
        x.pr[t,] <- colMeans(x_pr)
        
        # [2] (Update step)
        lik <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x_pr[p,], sigma=(var.eps * diag(D)))) 
        log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
 
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # compute filtered estimator to update (a posteriori) state estimate
        x_up <- m.sir(x_pr, lik, u.sim[t,])
        x.up.particles[[t]] <- x_up
        x.up[t,] <- colMeans(x_up) 
    }
    
    return(list(x.pr=x.pr, x.up=x.up, x.pr.particles=x.pr.particles, x.up.particles=x.up.particles, loglik=loglik))
}

# ----------------------------------------------------------------------
# (Univariate) maximum likelihood estimator
# ----------------------------------------------------------------------
particle.mle <- function(y, P, verbose=TRUE) {
    T <- length(y)
    
    # draw noise and particles
    eta.sim <- matrix(rnorm(T*P, mean=0, sd=1), nrow=T, ncol=P) 
    u.sim   <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
    for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
    
    # set optimization parameters
    lb <- 0.1
    ub <- 5
    theta_start <- 2
    obj <- function(theta){ return( -particle.filter(y, cov.eta=theta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=FALSE)$loglik ) } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- particle.filter(y, cov.eta=theta_mle, eta.sim=eta.sim, u.sim=u.sim)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}

# ----------------------------------------------------------------------
# (Multivariate) maximum likelihood estimator
# ----------------------------------------------------------------------
# note: does not work properly for multivariate case
m.particle.mle <- function(y, D, P, verbose=TRUE) {
    T <- nrow(y)
    
    # draw noise and particles
    eta.sim <- list()
    for (t in 1:T) {
        eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D)) 
    }
    u.sim <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
    for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
    
    # set optimization parameters
    lb <- c(rep(0.1,D), -1)
    ub <- c(rep(5,  D),  1)
    theta_start <- c(rep(1,D), 0)
    obj <- function(theta){ return( -m.particle.filter(y, cov.eta=construct.cov(c(theta[1:D]), theta[D+1]), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik ) } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- m.particle.filter(y, cov.eta=construct.cov(theta_mle[1:D], theta_mle[D+1]), eta.sim=eta.sim, u.sim=u.sim)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
