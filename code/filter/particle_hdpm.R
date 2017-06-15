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
    #x.pr.particles <- x.up.particles <- matrix(nrow = N*M, ncol = P)
    x.pr.particles <- x.up.particles <- list()
    loglik <- 0 
    
    # initial state estimator
    x_up.D <- x_up.init
    x_up.I <- x_up.init
    
    # estimate states of state space model
    for (n in 1:N) {    
        for (m in 1:M) {
            t <- (n-1)*M + m
  
            # [1] (Prediction step) 
            # one-step ahead (a priori) prediction          
            x_pr.D <- theta$D.phi0 + theta$D.phi1 * x_up.D + sqrt(theta$D.var) * noise.sim$D[n,]
            x_pr.P <- sin(pi*(m-1)/M) * theta$P.int
            x_pr.I <- theta$I.phi1 * x_up.I + sqrt(theta$I.var) * noise.sim$I[t,]
            x_pr   <- x_pr.D + x_pr.P + x_pr.I
        
            #x.pr.particles[t,] <- exp(x_pr)
            x.pr.particles[[t]] <- cbind(x_pr.D, rep(x_pr.P,P), x_pr.I)
            #x.pr[t] <- mean(exp(x_pr))
            x.pr[t] <- exp(mean(x_pr))
            
            
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
            #x.up.particles[t,] <- exp(x_up)
            x.up.particles[[t]] <- cbind(x_up.D, rep(x_up.P,P), x_up.I)
            #x.up[t] <- mean(exp(x_up)) 
            x.up[t] <- exp(mean(x_up)) 
        }
    }
    
    return(list(x.pr=x.pr, x.up=x.up, x.pr.particles=x.pr.particles, x.up.particles=x.up.particles, loglik=loglik))
}

# ----------------------------------------------------------------------
# Maximum likelihood estimator
# ----------------------------------------------------------------------
particle.mle.hdpm <- function(y, P, verbose=TRUE) {
    N <- nrow(y)
    M <- ncol(y)
    
    # draw standard normal transition noise 
    noise.sim <- list()
    noise.sim$D <- matrix(rnorm(N*P, mean=0, sd=1), nrow=N, ncol=P)
    noise.sim$I <- matrix(rnorm(N*M*P, mean=0, sd=1), nrow=N*M, ncol=P)
    
    # draw particles from standard uniform 
    u.sim   <- matrix(runif(N*M*P, min=0, max=1), nrow=N*M, ncol=P)   
    for (nm in c(1:(N*M))) {u.sim[nm,] <- sort( u.sim[nm,] )}
    
    # set optimization parameters
    # (D.phi0, D.phi1, I.phi1, P.int, D.var, I.var)
    lb <- c(0, 0, 0, 0, 0.1, 0.1)
    ub <- c(1, 1, 1, 1, 5, 5)
    theta_start <- c(0, 0, 0, 0, 1, 1)
    #theta_start <- list(D.phi0=0, D.phi1=0, I.phi1=0, P.int=0, D.var=1, I.var=1)
    obj <- function(theta) { 
        theta.l <- list(D.phi0=theta[1], D.phi1=theta[2], I.phi1=theta[3], P.int=theta[4], D.var=theta[5], I.var=theta[6])
        return( -particle.filter.hdpm(y, theta=theta.l, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik ) 
    } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    theta_mle.l <- list(D.phi0=theta_mle[1], D.phi1=theta_mle[2], I.phi1=theta_mle[3], P.int=theta_mle[4], D.var=theta_mle[5], I.var=theta_mle[6])
    loglik <- particle.filter.hdpm(y, theta=theta_mle.l, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle.l))
}
