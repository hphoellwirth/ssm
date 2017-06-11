# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Run space state methods on hierarchical dynamic Poisson model
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
source("models/hierarchDynPoisson.R")

# load filters
source("filter/particle_hdpm.R")
source("filter/aux.R")

# load utilities
source("util/plots.R")


# ----------------------------------------------------------------------
# Simulate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------

# generate hierarchical dynamic Poisson data
set.seed(1012)
N <- 5; M <- 20
theta <- list(D.phi0=0.7, D.phi1=0.6, I.phi1=0.3, P.int=0.8, D.var=0.6, I.var=0.2)
hdpm.data <- gen.hdpm.data(N=N, M=M, theta, a1=0, P1=1)

# plot data and parameter (components)
if(save.plots) png("../images/dyn-poisson.png", width=1000, height=500, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(as.vector(t(as.matrix(hdpm.data$x))), type='p', pch=19, col="red", xlab="time", ylab="counts / parameter (components)")
lines(as.vector(t(as.matrix(hdpm.data$lambda))), type='l', col="black")
lines(hdpm.data$D, type='l', col="orange")
lines(hdpm.data$P, type='l', col="green")
lines(hdpm.data$I, type='l', col="blue")
if(save.plots) dev.off()

# plot log parameter (components)
if(save.plots) png("../images/dyn-poisson-log-param.png", width=1000, height=500, pointsize=14)
plot(log(as.vector(t(as.matrix(hdpm.data$lambda)))), type='l', col="black", xlab="time", ylab="log parameter (components)")
lines(log(hdpm.data$D), type='l', col="orange")
lines(log(hdpm.data$P), type='l', col="green")
lines(log(hdpm.data$I), type='l', col="blue")
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Use filters with true parameters to estimate model states
# ----------------------------------------------------------------------
set.seed(1012)
P <- 200 
hdpm.particle.filter <- particle.filter.hdpm(hdpm.data$x, theta, P=P, x_up.init=rep(0,P))

# plot observations, states, and estimates
if(save.plots) png("../images/hdpm_est.png", width=1000, height=600, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(as.vector(t(as.matrix(hdpm.data$x))), type='p', pch=19, col="red", xlab="time", ylab="counts & parameter")
lines(as.vector(t(as.matrix(hdpm.data$lambda))), type='l', col="blue")
lines(hdpm.particle.filter$x.pr, col="orange")
#legend(10,15, c('observation','state','kalman','particle'), cex=1.0, lty=rep(1,4), lwd=rep(2.5,4), col=c('red','blue','green','orange'))
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Produce log-likelihood plots to validate and compare filters
# ----------------------------------------------------------------------
# draw noise and particles for comparison of particle filter
noise.sim <- list()
noise.sim$D <- matrix(rnorm(N*P, mean=0, sd=1), nrow=N, ncol=P)
noise.sim$I <- matrix(rnorm(N*M*P, mean=0, sd=1), nrow=N*M, ncol=P)

u.sim <- matrix(runif(N*M*P, min=0, max=1), nrow=N*M, ncol=P)   
for (nm in c(1:(N*M))) {u.sim[nm,] <- sort( u.sim[nm,] )}

# compute log-likelihood for wide range of different D.phi0 values
D.phi0 <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.D.phi0 <- rep(0,length(D.phi0))
for (i in 1:length(D.phi0)) {
    cat('.')
    theta.i$D.phi0 <- D.phi0[i]
    hdpm.particle.D.phi0[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

# compute log-likelihood for wide range of different D.phi1 values
D.phi1 <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.D.phi1 <- rep(0,length(D.phi1))
for (i in 1:length(D.phi1)) {
    cat('.')
    theta.i$D.phi1 <- D.phi1[i]
    hdpm.particle.D.phi1[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

# compute log-likelihood for wide range of different I.phi1 values
I.phi1 <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.I.phi1 <- rep(0,length(I.phi1))
for (i in 1:length(I.phi1)) {
    cat('.')
    theta.i$I.phi1 <- I.phi1[i]
    hdpm.particle.I.phi1[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

# compute log-likelihood for wide range of different I.phi1 values
P.int <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.P.int <- rep(0,length(P.int))
for (i in 1:length(P.int)) {
    cat('.')
    theta.i$P.int <- P.int[i]
    hdpm.particle.P.int[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

# compute log-likelihood for wide range of different D.var values
D.var <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.D.var <- rep(0,length(D.var))
for (i in 1:length(D.var)) {
    cat('.')
    theta.i$D.var <- D.var[i]
    hdpm.particle.D.var[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

# compute log-likelihood for wide range of different I.var values
I.var <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.I.var <- rep(0,length(I.var))
for (i in 1:length(D.var)) {
    cat('.')
    theta.i$I.var <- I.var[i]
    hdpm.particle.I.var[i] <- particle.filter.hdpm(hdpm.data$x, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

if(save.plots) png("../images/hdpm-loglik-1.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(D.phi0, hdpm.particle.D.phi0, theta$D.phi0, 'orange', 'D.phi0 with particle filter')
plot.loglik(D.phi1, hdpm.particle.D.phi1, theta$D.phi1, 'orange', 'D.phi1 with particle filter')
plot.loglik(I.phi1, hdpm.particle.I.phi1, theta$I.phi1, 'orange', 'I.phi1 with particle filter')
if(save.plots) dev.off()

if(save.plots) png("../images/hdpm-loglik-2.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(P.int, hdpm.particle.P.int, theta$P.int, 'orange', 'P.int with particle filter')
plot.loglik(D.var, hdpm.particle.D.var, theta$D.var, 'orange', 'D.var with particle filter')
plot.loglik(I.var, hdpm.particle.I.var, theta$I.var, 'orange', 'I.var with particle filter')
if(save.plots) dev.off()




