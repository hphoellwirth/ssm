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
# (Univariate) Auxiliary filter
# ----------------------------------------------------------------------
aux.filter <- function(y, x.pr, x.up, var.eps=1, var.eps.aux=1, cov.eta=1, cov.eta.aux=1) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- ncol(x.pr)
    
    # initialize importance weights and loglikelihood
    is.pr <- matrix(nrow = T, ncol = P)
    is.up <- matrix(nrow = T, ncol = P)     
    loglik <- 0 

    # initial computation of predictive importance weights
    p     <- dnorm( x.pr[1,], mean=0, sd=sqrt(cov.eta))
    p.aux <- dnorm( x.pr[1,], mean=0, sd=sqrt(cov.eta.aux))
    is.pr[1,] <- p / p.aux   
    
    for (t in 1:T) {
        # estimate likelihood and auxiliary likelihood
        lik     <- dnorm( y[t,]*rep(1,P), mean=x.pr[t,], sd=sqrt(var.eps)) * is.pr[t,]
        lik.aux <- dnorm( y[t,]*rep(1,P), mean=x.pr[t,], sd=sqrt(var.eps.aux)) 
        
        log.mean.lik     <- tryCatch(log(mean(lik)), error=function(e)(-Inf))     
        log.mean.lik.aux <- tryCatch(log(mean(lik.aux)), error=function(e)(-Inf))
        
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # [1] update filtering importance weights
        m     <- dnorm( y[t,]*rep(1,P), mean=x.up[t,], sd=sqrt(var.eps))  
        m.aux <- dnorm( y[t,]*rep(1,P), mean=x.up[t,], sd=sqrt(var.eps.aux))
        is.up[t,] <- (m / m.aux) * (mean(lik.aux) / mean(lik)) * is.pr[t,][match(x.up[t,], x.pr[t,])] 
        
        # [2] update predictive importance weights
        if (t < T) { 
            p     <- dnorm( x.pr[(t+1),], mean=x.up[t,], sd=sqrt(cov.eta))
            p.aux <- dnorm( x.pr[(t+1),], mean=x.up[t,], sd=sqrt(cov.eta.aux))
            is.pr[(t+1),] <- p / p.aux * is.up[t,]   
        }
    }
    
    return(list(is.pr=is.pr, is.up=is.up, loglik=loglik))
}

# ----------------------------------------------------------------------
# (Multivariate) Auxiliary filter
# ----------------------------------------------------------------------
m.aux.filter <- function(y, D=ncol(data.frame(y)), x.pr, x.up, var.eps=1, var.eps.aux=1, cov.eta=diag(D), cov.eta.aux=diag(D)) {
    y <- data.frame(y)
    T <- nrow(y)
    P <- nrow(x.pr[[1]])
    
    cov.eps     <- var.eps * diag(D)
    cov.eps.aux <- var.eps.aux * diag(D)
    
    # initialize importance weights and loglikelihood
    is.pr <- matrix(nrow = T, ncol = P)
    is.up <- matrix(nrow = T, ncol = P)     
    loglik <- 0 
    
    # initial computation of predictive importance weights
    p     <- sapply(1:P, function(p) dmvnorm(x.pr[[1]][p,], mean=rep(0,D), sigma=cov.eta)) 
    p.aux <- sapply(1:P, function(p) dmvnorm(x.pr[[1]][p,], mean=rep(0,D), sigma=cov.eta.aux)) 
    is.pr[1,] <- p / p.aux   
    
    for (t in 1:T) {
        # estimate likelihood and auxiliary likelihood
        lik     <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x.pr[[t]][p,], sigma=cov.eps)) * is.pr[t,] 
        lik.aux <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x.pr[[t]][p,], sigma=cov.eps.aux)) 
        
        log.mean.lik     <- tryCatch(log(mean(lik)), error=function(e)(-Inf))     
        log.mean.lik.aux <- tryCatch(log(mean(lik.aux)), error=function(e)(-Inf))
               
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # [1] update filtering importance weights
        m     <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x.up[[t]][p,], sigma=cov.eps)) 
        m.aux <- sapply(1:P, function(p) dmvnorm(y[t,], mean=x.up[[t]][p,], sigma=cov.eps.aux)) 
        is.up[t,] <- (m / m.aux) * (mean(lik.aux) / mean(lik)) * is.pr[t,][match.rows(x.up[[t]], x.pr[[t]])] 
                
        # [2] update predictive importance weights
        if (t < T) { 
            p     <- sapply(1:P, function(p) dmvnorm(x.pr[[t+1]][p,], mean=x.up[[t]][p,], sigma=cov.eta)) 
            p.aux <- sapply(1:P, function(p) dmvnorm(x.pr[[t+1]][p,], mean=x.up[[t]][p,], sigma=cov.eta.aux)) 
            is.pr[(t+1),] <- p / p.aux * is.up[t,]   
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
    theta_aux <- c(rep(1.5,D),0.2)
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
    theta_start <- c(rep(0.5,D), 0.5)
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
