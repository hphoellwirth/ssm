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
    P <- nrow(u.sim)
    
    # initialize importance weights and loglikelihood
    is.pr <- matrix(nrow = P, ncol = T)
    is.up <- matrix(nrow = P, ncol = T)     
    loglik <- 0 
    
    # initial computation of predictive importance weights
    p     <- dnorm( x.pr[,1], mean=0, sd=sqrt(cov.eta))
    p.aux <- dnorm( x.pr[,1], mean=0, sd=sqrt(cov.eta.aux))
    is.pr[,1] <- p / p.aux   
    
    # 
    for (t in 1:T) {
        # estimate likelihood and auxiliary likelihood
        lik     <- dnorm( y[t,]*rep(1,P), mean=x.pr[,t], sd=sqrt(var.eps)) * is.pr[,t]
        lik.aux <- dnorm( y[t,]*rep(1,P), mean=x.pr[,t], sd=sqrt(var.eps.aux)) 
        
        log.mean.lik     <- tryCatch(log(mean(lik)), error=function(e)(-Inf))     
        log.mean.lik.aux <- tryCatch(log(mean(lik.aux)), error=function(e)(-Inf))
        
        # update likelihood
        if (log.mean.lik == -Inf) {
            loglik <- -Inf
            break
        }
        loglik <- loglik + log.mean.lik
        
        # [1] update filtering importance weights
        m     <- dnorm( y[t,]*rep(1,P), mean=x.up[,t], sd=sqrt(var.eps))  
        m.aux <- dnorm( y[t,]*rep(1,P), mean=x.up[,t], sd=sqrt(var.eps.aux))
        is.up[,t] <- (m / m.aux) * (mean(lik.aux) / mean(lik)) * is.pr[,t][match(x.up[,t], x.pr[,t])] 
        
        # [2] update predictive importance weights
        if (t < T) { 
            p     <- dnorm( x.pr[,(t+1)], mean=x.up[,t], sd=sqrt(cov.eta))
            p.aux <- dnorm( x.pr[,(t+1)], mean=x.up[,t], sd=sqrt(cov.eta.aux))
            is.pr[,(t+1)] <- p / p.aux * is.up[,t]   
        }
    }
    
    return(list(loglik=loglik))
}

# ----------------------------------------------------------------------
# (Univariate) maximum likelihood estimator
# ----------------------------------------------------------------------
aux.mle <- function(y, P) {
    T <- length(y)
    
    # draw noise and particles
    theta_aux = 1
    eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=P, ncol=T) 
    u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
    for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
    pf <- particle.filter(y, cov.eta=theta_aux, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))
    
    # set optimization parameters
    lb <- 0.1
    ub <- 5
    theta_start <- 1
    obj <- function(theta){ return( -aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=theta, cov.eta.aux=theta_aux)$loglik ) } 
    
    # run box-constrained optimization
    print('estimating model parameters...') 
    param <- nlminb( theta_start, obj, lower=lb, upper=ub )
    theta_mle <- param$par
    print('... done!') 
    
    # compute log-liklihood of MLE parameters
    loglik <- aux.filter(y, x.pr=pf$x.pr.particles, x.up=pf$x.up.particles, cov.eta=theta_mle, cov.eta.aux=theta_aux)$loglik
    
    return(list(loglik = loglik, theta_mle = theta_mle))
}
