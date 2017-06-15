# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Auxiliary filter (for hierarchical dynamic Poisson model)
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
# Maximum likelihood estimator
# ----------------------------------------------------------------------
aux.mle.hdpm <- function(y, P, verbose=TRUE) {
    N <- nrow(y)
    M <- ncol(y)
    
    # draw standard normal transition noise 
    noise.sim <- list()
    noise.sim$D <- matrix(rnorm(N*P, mean=0, sd=1), nrow=N, ncol=P)
    noise.sim$I <- matrix(rnorm(N*M*P, mean=0, sd=1), nrow=N*M, ncol=P)
    
    # draw particles from standard uniform 
    u.sim   <- matrix(runif(N*M*P, min=0, max=1), nrow=N*M, ncol=P)   
    for (nm in c(1:(N*M))) {u.sim[nm,] <- sort( u.sim[nm,] )}    
    theta_aux <- 1
    eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=T, ncol=P) 

    theta_aux <- list(D.phi0=0.5, D.phi1=0.5, I.phi1=0.5, P.int=1, D.var=1, I.var=1)
    pf <- particle.filter.hdpm(y, theta_aux, P=P, x_up.init=rep(0,P))
    
    # set optimization parameters
    # (D.phi0, D.phi1, I.phi1, P.int, D.var, I.var)
    lb <- c(0, 0, 0, 0, 0.1, 0.1)
    ub <- c(1, 1, 1, 1, 5, 5)
    theta_start <- c(0, 0, 0, 0, 1, 1)
    obj <- function(theta) { 
        theta.l <- list(D.phi0=theta[1], D.phi1=theta[2], I.phi1=theta[3], P.int=theta[4], D.var=theta[5], I.var=theta[6])
        return( -aux.filter.hdpm(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, theta=theta.l, theta.aux=theta_aux)$loglik ) 
    } 
    
    # run box-constrained optimization
    if (verbose) print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    if (verbose) print('... done!') 
    
    # compute log-liklihood of MLE parameters
    theta_mle.l <- list(D.phi0=theta_mle[1], D.phi1=theta_mle[2], I.phi1=theta_mle[3], P.int=theta_mle[4], D.var=theta_mle[5], I.var=theta_mle[6])
    loglik <- aux.filter.hdpm(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, theta=theta_mle.l, theta.aux=theta_aux)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle.l))
}


