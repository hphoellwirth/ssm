# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Run space state methods on multivariate local level model
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# house cleaning
rm(list = ls())
par(mfrow=c(1,1))
save.plots <- FALSE

# load libraries
# NONE

# if interactive, during the development, set to TRUE
interactive <- FALSE
if (interactive) {
    setwd("/Users/Hans-Peter/Documents/Masters/14D000/code")
} 

# load model
source("models/localLevel.R")

# load filters
source("filter/kalman.R")
source("filter/particle.R")
source("filter/aux.R")

# load utilities
source("util/plots.R")
source("util/misc.R")


# ----------------------------------------------------------------------
# Generate trivariate local level data
# ----------------------------------------------------------------------
set.seed(1000)
D <- 3
T <- 50
cov.eta.var <- c(4.2,2.8,0.9)
cov.eta.rho <- 0.7
mllm.data <- gen.multi.llm.data(n=T, d=D, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))


# ----------------------------------------------------------------------
# Use filters with true parameters to estimate model states
# ----------------------------------------------------------------------

# use Kalman filter to estimate model states
mllm.kalman.filter <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use particle filter to estimate model states
P <- 200
eta.sim <- list()
for (t in 1:T) {eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D))}
u.sim <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
mllm.particle.filter <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))

# plot observations, states, and estimates
if(save.plots) png("../images/trivariate-local-level.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(2,2,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red")
    lines(mllm.data$x[,d], col="blue")
    lines(mllm.kalman.filter$a[,d], col="green")
    lines(mllm.particle.filter$x.pr[,d], col="orange")
}
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Evaluate the correctness of the multivariate auxiliary filter
# ----------------------------------------------------------------------

# use auxiliary filter to compute true log-likelihood
mllm.aux.filter <- m.aux.filter(mllm.data$y, x.pr=mllm.particle.filter$x.pr.particles, x.up=mllm.particle.filter$x.up.particles, cov.eta=construct.cov(cov.eta.var, cov.eta.rho), cov.eta.aux=diag(D))

# plot importance weights over time
if(save.plots) png("../images/ullm_aux_weights.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.weights(mllm.aux.filter$is.up, xlab=paste0('time T (with P=',P,' particles)'), ylab='filtering weights')
plot.weights(mllm.aux.filter$is.pr, xlab=paste0('time T (with P=',P,' particles)'), ylab='predictive weights')
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Produce log-likelihood plots to compare filters
# ----------------------------------------------------------------------

# draw standard normal transition noise and particles for comparison of particle filter
eta.sim <- list()
for (t in 1:T) {eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D))}
u.sim <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}

# plot log-likelihood for different rho values
rho <- seq(-1.0,1,0.1)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(rho))
for (i in 1:length(rho)) {
    cat('.')
    cov.eta        <- construct.cov(cov.eta.var, rho[i])
    ll.kalman[i]   <- kalman.filter(mllm.data$y, cov.eta=cov.eta)$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=cov.eta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- m.aux.filter(mllm.data$y, x.pr=mllm.particle.filter$x.pr.particles, x.up=mllm.particle.filter$x.up.particles, cov.eta=cov.eta, cov.eta.aux=construct.cov(c(1,1,1),0))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-rho.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(rho, ll.kalman, cov.eta.rho, 'green', 'rho with Kalman filter')
plot.loglik(rho, ll.particle, cov.eta.rho, 'orange', 'rho with particle filter')
plot.loglik(rho, ll.aux, cov.eta.rho, 'magenta', 'rho with auxiliary filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var1 values
var1 <- seq(2,5.6,0.2)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(var1))
for (i in 1:length(var1)) {
    cat('.')
    cov.eta        <- construct.cov(c(var1[i],cov.eta.var[2],cov.eta.var[3]), cov.eta.rho)
    ll.kalman[i]   <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(var1[i],2.8,0.9), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=cov.eta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- m.aux.filter(mllm.data$y, x.pr=mllm.particle.filter$x.pr.particles, x.up=mllm.particle.filter$x.up.particles, cov.eta=cov.eta, cov.eta.aux=construct.cov(c(1,1,1),0))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var1.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(var1, ll.kalman, cov.eta.var[1], 'green', 'eta.var1 with Kalman filter')
plot.loglik(var1, ll.particle, cov.eta.var[1], 'orange', 'eta.var1 with particle filter')
plot.loglik(var1, ll.aux, cov.eta.var[1], 'magenta', 'eta.var1 with auxiliary filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var2 values
var2 <- seq(1.6,4.6,0.2)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(var2))
for (i in 1:length(var2)) {
    cat('.')
    cov.eta        <- construct.cov(c(cov.eta.var[1],var2[i],cov.eta.var[3]), cov.eta.rho)
    ll.kalman[i]   <- kalman.filter(mllm.data$y, cov.eta=cov.eta)$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=cov.eta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- m.aux.filter(mllm.data$y, x.pr=mllm.particle.filter$x.pr.particles, x.up=mllm.particle.filter$x.up.particles, cov.eta=cov.eta, cov.eta.aux=construct.cov(c(1,1,1),0))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var2.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(var2, ll.kalman, cov.eta.var[2], 'green', 'eta.var2 with Kalman filter')
plot.loglik(var2, ll.particle, cov.eta.var[2], 'orange', 'eta.var2 with particle filter')
plot.loglik(var2, ll.aux, cov.eta.var[2], 'magenta', 'eta.var2 with auxiliary filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var3 values
var3 <- seq(0.1,3.9,0.2)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(var3))
for (i in 1:length(var3)) {
    cat('.')
    cov.eta        <- construct.cov(c(cov.eta.var[1],cov.eta.var[2],var3[i]), cov.eta.rho)
    ll.kalman[i]   <- kalman.filter(mllm.data$y, cov.eta=cov.eta)$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=cov.eta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- m.aux.filter(mllm.data$y, x.pr=mllm.particle.filter$x.pr.particles, x.up=mllm.particle.filter$x.up.particles, cov.eta=cov.eta, cov.eta.aux=construct.cov(c(1,1,1),0))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var3.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(var3, ll.kalman, cov.eta.var[3], 'green', 'eta.var3 with Kalman filter')
plot.loglik(var3, ll.particle, cov.eta.var[3], 'orange', 'eta.var3 with particle filter')
plot.loglik(var3, ll.aux, cov.eta.var[3], 'magenta', 'eta.var3 with auxiliary filter')
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Compute MLE with different filters
# ----------------------------------------------------------------------

# estimate model parameters, using Kalman filter
mllm.kalman.mle <- kalman.mle(mllm.data$y, D=D)
print(paste('     True parameters:', round(cov.eta.var[1],2), round(cov.eta.var[2],2), round(cov.eta.var[3],2), round(cov.eta.rho,2)))
print(paste('Estimated parameters:', round(mllm.kalman.mle$theta_mle[1],2), round(mllm.kalman.mle$theta_mle[2],2), round(mllm.kalman.mle$theta_mle[3],2), round(mllm.kalman.mle$theta_mle[4],2)))
print(paste('     True log-likelihood:', round(mllm.kalman.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(mllm.kalman.mle$loglik,3)))

# estimate model parameters, using particle filter
mllm.particle.mle <- m.particle.mle(mllm.data$y, D=D, P=200)
print(paste('     True parameters:', round(cov.eta.var[1],2), round(cov.eta.var[2],2), round(cov.eta.var[3],2), round(cov.eta.rho,2)))
print(paste('Estimated parameters:', round(mllm.particle.mle$theta_mle[1],2), round(mllm.particle.mle$theta_mle[2],2), round(mllm.particle.mle$theta_mle[3],2), round(mllm.particle.mle$theta_mle[4],2)))
print(paste('     True log-likelihood:', round(mllm.particle.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(mllm.particle.mle$loglik,3)))

# plot Kalman filter with true and estimated parameters
if(save.plots) png("../images/trivariate-local-level-est-kalman.png", width=1000, height=750, pointsize=20)
mllm.kalman.est <- kalman.filter(mllm.data$y, cov.eta=construct.cov(mllm.kalman.mle$theta_mle[1:3], mllm.kalman.mle$theta_mle[4]))
par(mfrow=c(3,1), mar=c(2,2,1,1))
for (d in 1:3) {
    plot(mllm.data$x[,d], type='l', col="blue") 
    lines(mllm.kalman.filter$a[,d], col="green")
    lines(mllm.kalman.est$a[,d], col="purple")
}
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Compare filters with true and estimated (MLE) parameters
# ----------------------------------------------------------------------

# plot particle filter with true and estimated parameters
if(save.plots) png("../images/trivariate-local-level-est-particle.png", width=1000, height=750, pointsize=20)
mllm.particle.est <- particle.filter(mllm.data$y, cov.eta=construct.cov(mllm.particle.mle$theta_mle[1:3], mllm.particle.mle$theta_mle[4]), eta.sim=eta.sim, u.sim=u.sim)
par(mfrow=c(3,1), mar=c(2,2,1,1))
for (d in 1:3) {
    plot(mllm.data$x[,d], type='l', col="blue") 
    lines(mllm.particle.filter$x.pr[,d], col="green")
    lines(mllm.particle.est$x.pr[,d], col="purple")
}
if(save.plots) dev.off()

# and compare MSE of true and estimated parameters
mse <- function(target, pred) {
    return(sum((target - pred)**2))
}
print(paste('     True MSE:', mse(mllm.data$x, mllm.filter$a)))
print(paste('Estimated MSE:', mse(mllm.data$x, mllm.est$a)))





