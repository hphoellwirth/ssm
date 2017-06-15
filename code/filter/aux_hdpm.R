# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Auxiliary filter (for multivariate local level model)
#
# (Author) Hans-Peter Höllwirth
# (Date)   30.05.2017

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
library(matrixcalc)
library(mvtnorm)
library(nloptr)

# ----------------------------------------------------------------------
# Auxiliary filter
# ----------------------------------------------------------------------
aux.filter.hdpm <- function(y, x.pr, x.up, theta, theta.aux) {
    y <- data.frame(y)
    N <- nrow(y)
    M <- ncol(y)
    P <- nrow(x.pr[[1]])
    
    # initialize importance weights and loglikelihood
    is.pr <- matrix(nrow = N*M, ncol = P)
    is.up <- matrix(nrow = N*M, ncol = P)     
    loglik <- 0 

    # initial computation of predictive importance weights
    x_pr.D <- theta$D.phi0 + theta$D.phi1 * 0
    x_pr.P <- sin(pi*0/M) * theta$P.int
    x_pr.I <- theta$I.phi1 * 0
    x_pr   <- rep(x_pr.D + x_pr.P + x_pr.I, P)
    
    x_pr.aux.D <- theta.aux$D.phi0 + theta.aux$D.phi1 * 0
    x_pr.aux.P <- sin(pi*0/M) * theta.aux$P.int
    x_pr.aux.I <- theta.aux$I.phi1 * 0
    x_pr.aux   <- rep(x_pr.aux.D + x_pr.aux.P + x_pr.aux.I, P)
    
    p     <- dnorm( rowSums(x.pr[[1]]), mean=x_pr, sd=sqrt(theta$D.var + theta$I.var)) 
    p.aux <- dnorm( rowSums(x.pr[[1]]), mean=x_pr.aux, sd=sqrt(theta.aux$D.var + theta.aux$I.var))    
    is.pr[1,] <- p / p.aux   
    
    for (n in 1:N) {
        for (m in 1:M) {
            t <- (n-1)*M + m
            x_pr.sum <- rowSums(x.pr[[t]])
            x_up.sum <- rowSums(x.up[[t]])
            
            # estimate likelihood and auxiliary likelihood
            lik     <- dpois( y[n,m]*rep(1,P), lambda=exp(x_pr.sum)) * is.pr[t,]
            lik.aux <- dpois( y[n,m]*rep(1,P), lambda=exp(x_pr.sum)) 
            
            log.mean.lik     <- tryCatch(log(mean(lik)), error=function(e)(-Inf))     
            log.mean.lik.aux <- tryCatch(log(mean(lik.aux)), error=function(e)(-Inf))
            
            # update likelihood
            if (log.mean.lik == -Inf) {
                loglik <- -Inf
                break
            }
            loglik <- loglik + log.mean.lik
            
            # [1] update filtering importance weights
            em     <- dpois( y[n,m]*rep(1,P), lambda=exp(x_up.sum))
            em.aux <- dpois( y[n,m]*rep(1,P), lambda=exp(x_up.sum))
            is.up[t,] <- (em / em.aux) * (mean(lik.aux) / mean(lik)) * is.pr[t,][match(x_up.sum, x_pr.sum)] 
            
            # [2] update predictive importance weights
            if (t < (N*M)) { 
                x_pr.D <- theta$D.phi0 + theta$D.phi1 * x.up[[t]][,1]
                x_pr.P <- sin(pi*(m-1)/M) * theta$P.int
                x_pr.I <- theta$I.phi1 * x.up[[t]][,3]
                x_pr   <- x_pr.D + x_pr.P + x_pr.I
                
                x_pr.aux.D <- theta.aux$D.phi0 + theta.aux$D.phi1 * x.up[[t]][,1]
                x_pr.aux.P <- sin(pi*(m-1)/M) * theta.aux$P.int
                x_pr.aux.I <- theta.aux$I.phi1 * x.up[[t]][,3]
                x_pr.aux   <- x_pr.aux.D + x_pr.aux.P + x_pr.aux.I                
                
                p     <- dnorm( rowSums(x.pr[[t+1]]), mean=x_pr, sd=sqrt(theta$D.var + theta$I.var)) 
                p.aux <- dnorm( rowSums(x.pr[[t+1]]), mean=x_pr.aux, sd=sqrt(theta.aux$D.var + theta.aux$I.var))
                is.pr[(t+1),] <- p / p.aux * is.up[t,]   
            }
        }
    }
    
    return(list(is.pr=is.pr, is.up=is.up, loglik=loglik))
}

# ----------------------------------------------------------------------
# (Univariate) maximum likelihood estimator
# ----------------------------------------------------------------------
aux.mle <- function(y, P, verbose=TRUE) {
    T <- length(y)
    
    # draw noise and particles
    theta_aux <- 1
    eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=T, ncol=P) 
    u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=T, ncol=P)   
    for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
    pf <- particle.filter(y, cov.eta=theta_aux, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=FALSE)
    
    # set optimization parameters
    lb <- 0.1
    ub <- 5
    theta_start <- 1
    obj <- function(theta){ return( -aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=theta, cov.eta.aux=theta_aux)$loglik ) } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=theta_mle, cov.eta.aux=theta_aux)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}

# ----------------------------------------------------------------------
# (Multivariate) maximum likelihood estimator
# ----------------------------------------------------------------------
m.aux.mle <- function(y, D, P, verbose=TRUE) {
    T <- nrow(y)
    
    # draw noise and particles
    theta_aux <- c(rep(1,D),0)
    eta.sim <- list()
    for (t in 1:T) {
        eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D)) 
    }    
    u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=T, ncol=P)   
    for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
    pf <- m.particle.filter(y, cov.eta=construct.cov(c(theta_aux[1:D]), theta_aux[D+1]), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))
    
    # set optimization parameters
    lb <- c(rep(0.1,D), -1)
    ub <- c(rep(5,  D),  1)
    theta_start <- c(rep(1,D), 0)
    obj <- function(theta){ return( -m.aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=construct.cov(c(theta[1:D]), theta[D+1]), cov.eta.aux=construct.cov(c(theta_aux[1:D]), theta_aux[D+1]))$loglik ) } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- m.aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=construct.cov(c(theta_mle[1:D]), theta_mle[D+1]), cov.eta.aux=construct.cov(c(theta_aux[1:D]), theta_aux[D+1]))$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
