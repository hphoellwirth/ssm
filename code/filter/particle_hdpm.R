# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Particle filter (for hierarchical dynamic Poisson model)
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017

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
# Particle filter for hierarchical dynamic Poisson model
# ----------------------------------------------------------------------
particle.filter.hdpm <- function(y, theta, noise.sim=NA, u.sim=NA, P=ncol(u.sim), a1=0, P1=1, x_up.init=rnorm(P, mean=a1, sd=P1)) {
    y <- data.frame(y)
    N <- nrow(y)
    M <- ncol(y)
    
    # draw standard normal transition noise 
    if (is.na(noise.sim)[1][1]) {
        noise.sim <- list()
        noise.sim$D <- matrix(rnorm(N*P, mean=0, sd=1), nrow=N, ncol=P)
        noise.sim$I <- matrix(rnorm(N*M*P, mean=0, sd=1), nrow=N*M, ncol=P)
    }
    
    # draw particles from standard uniform 
    if (is.na(u.sim)[1][1]) {
        u.sim   <- matrix(runif(N*M*P, min=0, max=1), nrow=N*M, ncol=P)   
        for (nm in c(1:(N*M))) {u.sim[nm,] <- sort( u.sim[nm,] )}
    }
    
    # initialize series and loglikelihood
    x.pr <- x.up <- rep(0,N*M)
    x.pr.particles <- x.up.particles <- matrix(nrow = N*M, ncol = P)
    loglik <- 0 
    
    # initial state estimator
    x_up.D <- x_up.init
    x_up.I <- x_up.init
    
    # estimate states of state space model
    for (n in 1:N) {    
        # [1] (Prediction step) 
        # one-step ahead (a priori) prediction
        x_pr.D <- theta$D.phi0 + theta$D.phi1 * x_up.D + sqrt(theta$D.var) * noise.sim$D[n,]
        
        for (m in 1:M) {
            t <- (n-1)*M + m
            
            x_pr.P <- sin(pi*(m-1)/M) * theta$P.int
            x_pr.I <- theta$I.phi1 * x_up.I + sqrt(theta$I.var) * noise.sim$I[t,]
            x_pr   <- x_pr.D + x_pr.P + x_pr.I
        
            x.pr.particles[t,] <- exp(x_pr)
            x.pr[t] <- mean(exp(x_pr))
        
            # [2] (Update step)
            # compute importance weights
            lik <- dpois( y[n,m]*rep(1,P), lambda=exp(x_pr))
            log.mean.lik <- tryCatch(log(mean(lik)), error=function(e)(-Inf))
            
            # update likelihood
            if (log.mean.lik == -Inf) {
                loglik <- -Inf
                break
            }
            loglik <- loglik + log.mean.lik
            
            # compute filtered estimator to update (a posteriori) state estimate
            x_up <- m.sir(matrix(c(x_pr.D, x_pr.I), nrow=P, ncol=2), lik, u.sim[t,])
            
            x_up.D <- x_up[,1]
            x_up.I <- x_up[,2]
            x_up.P <- sin(pi*(m-1)/M) * theta$P.int
            
            x_up <- x_up.D + x_up.P + x_up.I
            x.up.particles[t,] <- exp(x_up)
            x.up[t] <- mean(exp(x_up)) 
        }
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
