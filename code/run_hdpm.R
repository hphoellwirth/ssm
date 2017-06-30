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
library(matrixStats)
library(scales)

# if interactive, during the development, set to TRUE
interactive <- FALSE
if (interactive) {
    setwd("/Users/Hans-Peter/Documents/Masters/14D000/code")
} 

# load model
source("models/hierarchDynPoisson.R")

# load filters
source("filter/particle_hdpm.R")
source("filter/aux_hdpm.R")

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
if(save.plots) png("../images/hdpm-realization.png", width=600, height=350, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(as.vector(t(as.matrix(hdpm.data$y))), type='p', pch=19, col=alpha("red",0.7), xlab="time", ylab="counts & parameter")
lines(as.vector(t(as.matrix(hdpm.data$lambda))), type='l', col="black")
legend('topleft', legend=c('observation','state'), col=c('red','black'), cex=1.0, lty=c(0,1), lwd=c(0,2.5), pch=c(19,NA))
if(save.plots) dev.off()

# plot log parameter (components)
if(save.plots) png("../images/hdpm-log-param.png", width=600, height=350, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
ylim <- c(min(log(hdpm.data$I)), max(log(as.vector(t(as.matrix(hdpm.data$lambda))))))
plot(log(as.vector(t(as.matrix(hdpm.data$lambda)))), type='l', ylim=ylim, col="black", xlab="time", ylab="log parameter (components)")
lines(log(hdpm.data$D), type='l', col="red")
lines(log(hdpm.data$I), type='l', col="blue")
lines(log(hdpm.data$P), type='l', col="purple")
legend('topleft', legend=c('state','daily','intra-daily','periodic'), col=c('black','red','blue','purple'), cex=1.0, lty=1, lwd=2.5)
if(save.plots) dev.off()



# ----------------------------------------------------------------------
# Use particle filter with true parameters to estimate model states
# ----------------------------------------------------------------------
set.seed(2000)
P <- 200 
hdpm.particle.filter <- particle.filter.hdpm(hdpm.data$y, theta, P=P, x_up.init=rep(0,P))

# plot observations, states, estimates, and confidence intervals
states <- matrix(nrow=N*M, ncol=P)
for (t in 1:(N*M)) {states[t,] <- exp(rowSums(hdpm.particle.filter$x.up.particles[[t]]))}

if(save.plots) png("../images/hdpm-est.png", width=700, height=400, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(as.vector(t(as.matrix(hdpm.data$y))), type='p', pch=19, col=alpha("red",0.7), xlab="time", ylab="counts & parameter")
lines(as.vector(t(as.matrix(hdpm.data$lambda))), type='l', col="black")
lines(hdpm.particle.filter$x.up, col="orange", lwd=1.5)
lines(rowQuantiles(states, probs=c(0.95)), col=alpha('orange',1.0), lty=2)
lines(rowQuantiles(states, probs=c(0.05)), col=alpha('orange',1.0), lty=2)
legend('topleft', c('observation','state','estimate','90% CI'), cex=1.0, lty=c(0,1,1,2), lwd=rep(2.5,4), col=c('red','black','orange','orange'), pch=c(19,NA,NA,NA))
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Evaluate the correctness of the auxiliary filter
# ----------------------------------------------------------------------

# use auxiliary filter to compute true log-likelihood
theta.aux <- list(D.phi0=0.5, D.phi1=0.5, I.phi1=0.5, P.int=1, D.var=1, I.var=1)
hdpm.aux.filter <- aux.filter.hdpm(hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta, theta.aux=theta.aux)

# plot importance weights over time
if(save.plots) png("../images/hdpm_aux_weights.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.weights(hdpm.aux.filter$is.up, xlab=paste0('time T (with P=',P,' particles)'), ylab='filtering weights')
plot.weights(hdpm.aux.filter$is.pr, xlab=paste0('time T (with P=',P,' particles)'), ylab='predictive weights')
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Produce log-likelihood plots to validate and compare filters
# ----------------------------------------------------------------------
# draw noise and particles for comparison of particle filter
set.seed(2000)
noise.sim <- list()
noise.sim$D <- matrix(rnorm(N*P, mean=0, sd=1), nrow=N, ncol=P)
noise.sim$I <- matrix(rnorm(N*M*P, mean=0, sd=1), nrow=N*M, ncol=P)

u.sim <- matrix(runif(N*M*P, min=0, max=1), nrow=N*M, ncol=P)   
for (nm in c(1:(N*M))) {u.sim[nm,] <- sort( u.sim[nm,] )}

# compute log-likelihood for wide range of different D.phi0 values
D.phi0 <- seq(0.1,1.0,0.05)
theta.i <- theta
hdpm.particle.D.phi0 <- hdpm.aux.D.phi0 <- rep(0,length(D.phi0))
for (i in 1:length(D.phi0)) {
    cat('.')
    theta.i$D.phi0 <- D.phi0[i]
    hdpm.particle.D.phi0[i] <- particle.filter.hdpm( hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.D.phi0[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Dphi0.png", width=1000, height=500, pointsize=15)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(D.phi0, hdpm.particle.D.phi0, theta$D.phi0, 'orange', 'D.phi0 with particle filter')
plot.loglik(D.phi0, hdpm.aux.D.phi0, theta$D.phi0, 'magenta', 'D.phi0 with auxiliary filter')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different D.phi1 values
D.phi1 <- seq(0.1,1.0,0.05)
theta.i <- theta
hdpm.particle.D.phi1 <- hdpm.aux.D.phi1 <- rep(0,length(D.phi1))
for (i in 1:length(D.phi1)) {
    cat('.')
    theta.i$D.phi1 <- D.phi1[i]
    hdpm.particle.D.phi1[i] <- particle.filter.hdpm( hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.D.phi1[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Dphi1.png", width=500, height=300, pointsize=14)
par(mfrow=c(2,1), mar=c(2,4,1,1))
plot.loglik(D.phi1, hdpm.particle.D.phi1, theta$D.phi1, 'orange', 'D.phi1 with particle filter')
plot.loglik(D.phi1, hdpm.aux.D.phi1, theta$D.phi1, 'magenta', 'D.phi1 with auxiliary filter')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different I.phi1 values
I.phi1 <- seq(0.1,1.0,0.05)
theta.i <- theta
hdpm.particle.I.phi1 <- hdpm.aux.I.phi1 <- rep(0,length(I.phi1))
for (i in 1:length(I.phi1)) {
    cat('.')
    theta.i$I.phi1 <- I.phi1[i]
    hdpm.particle.I.phi1[i] <- particle.filter.hdpm(hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.I.phi1[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Iphi1.png", width=1000, height=500, pointsize=15)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(I.phi1, hdpm.particle.I.phi1, theta$I.phi1, 'orange', 'I.phi1 with particle filter')
plot.loglik(I.phi1, hdpm.aux.I.phi1, theta$I.phi1, 'magenta', 'I.phi1 with auxiliary filter')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different I.phi1 values
P.int <- seq(0.1,1.0,0.05)
theta.i <- theta
hdpm.particle.P.int <- hdpm.aux.P.int <- rep(0,length(P.int))
for (i in 1:length(P.int)) {
    cat('.')
    theta.i$P.int <- P.int[i]
    hdpm.particle.P.int[i] <- particle.filter.hdpm(hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.P.int[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Pint.png", width=1000, height=500, pointsize=15)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(P.int, hdpm.particle.P.int, theta$P.int, 'orange', 'P.int with particle filter')
plot.loglik(P.int, hdpm.aux.P.int, theta$P.int, 'magenta', 'P.int with auxiliary filter')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different D.var values
D.var <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.D.var <- hdpm.aux.D.var <- rep(0,length(D.var))
for (i in 1:length(D.var)) {
    cat('.')
    theta.i$D.var <- D.var[i]
    hdpm.particle.D.var[i] <- particle.filter.hdpm(hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.D.var[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Dvar.png", width=500, height=300, pointsize=14)
par(mfrow=c(2,1), mar=c(2,4,1,1))
plot.loglik(D.var, hdpm.particle.D.var, theta$D.var, 'orange', 'D.var with particle filter')
plot.loglik(D.var, hdpm.aux.D.var, theta$D.var, 'magenta', 'D.var with auxiliary filter')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different I.var values
I.var <- seq(0.1,2,0.1)
theta.i <- theta
hdpm.particle.I.var <- hdpm.aux.I.var <- rep(0,length(I.var))
for (i in 1:length(D.var)) {
    cat('.')
    theta.i$I.var <- I.var[i]
    hdpm.particle.I.var[i] <- particle.filter.hdpm(hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.I.var[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Ivar.png", width=1000, height=500, pointsize=15)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(I.var, hdpm.particle.I.var, theta$I.var, 'orange', 'I.var with particle filter')
plot.loglik(I.var, hdpm.aux.I.var, theta$I.var, 'magenta', 'I.var with auxiliary filter')
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Produce zoomed-in log-likelihood plots to validate and compare filters
# ----------------------------------------------------------------------

# compute log-likelihood for wide range of different D.phi1 values
D.phi1 <- seq(0.550,0.660,0.001)
theta.i <- theta
hdpm.particle.D.phi1 <- hdpm.aux.D.phi1 <- rep(0,length(D.phi1))
for (i in 1:length(D.phi1)) {
    cat('.')
    theta.i$D.phi1 <- D.phi1[i]
    hdpm.particle.D.phi1[i] <- particle.filter.hdpm( hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.D.phi1[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Dphi1-zoom.png", width=500, height=300, pointsize=14)
par(mfrow=c(2,1), mar=c(2,1,1,1))
plot.loglik.zoom(D.phi1, hdpm.particle.D.phi1, theta$D.phi1, 'orange', 'D.phi1 with particle filter', 'SIR')
plot.loglik.zoom(D.phi1, hdpm.aux.D.phi1, theta$D.phi1, 'magenta', 'D.phi1 with auxiliary filter', 'IS particle')
if(save.plots) dev.off()

# compute log-likelihood for wide range of different D.var values
D.var <- seq(0.550,0.650,0.001)
theta.i <- theta
hdpm.particle.D.var <- hdpm.aux.D.var <- rep(0,length(D.var))
for (i in 1:length(D.var)) {
    cat('.')
    theta.i$D.var <- D.var[i]
    hdpm.particle.D.var[i] <- particle.filter.hdpm(hdpm.data$y, theta.i, noise.sim=noise.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    hdpm.aux.D.var[i]      <- aux.filter.hdpm( hdpm.data$y, x.pr=hdpm.particle.filter$x.pr.particles, x.up=hdpm.particle.filter$x.up.particles, theta=theta.i, theta.aux=theta.aux)$loglik
}

if(save.plots) png("../images/hdpm-loglik-Dvar-zoom.png", width=500, height=300, pointsize=14)
par(mfrow=c(2,1), mar=c(2,1,1,1))
plot.loglik.zoom(D.var, hdpm.particle.D.var, theta$D.var, 'orange', 'D.var with particle filter','SIR')
plot.loglik.zoom(D.var, hdpm.aux.D.var, theta$D.var, 'magenta', 'D.var with auxiliary filter','IS particle')
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Compute MLE with different filters
# ----------------------------------------------------------------------
hdpm.mle.result <- data.frame(matrix(nrow=2, ncol=7), row.names=c('True','MLE'))
colnames(hdpm.mle.result) <- c('log-lik','D.phi0','D.phi1','I.phi1','P.int','D.var','I.var')

# estimate model parameters, using particle filter
hdpm.particle.mle <- particle.mle.hdpm(hdpm.data$y, P=200)
theta_mle <- hdpm.particle.mle$theta_mle

hdpm.particle.mle.result <- hdpm.mle.result
hdpm.particle.mle.result['True','log-lik'] <- round(hdpm.particle.filter$loglik,3)
hdpm.particle.mle.result['MLE', 'log-lik'] <- round(hdpm.particle.mle$loglik,3)
hdpm.particle.mle.result['True',2:7] <- round(c(theta$D.phi0, theta$D.phi1, theta$I.phi1, theta$P.int, theta$D.var, theta$I.var),2)
hdpm.particle.mle.result['MLE', 2:7] <- round(c(theta_mle$D.phi0, theta_mle$D.phi1, theta_mle$I.phi1, theta_mle$P.int, theta_mle$D.var, theta_mle$I.var),2)
hdpm.particle.mle.result

# estimate model parameters, using auxiliary filter
hdpm.aux.mle <- aux.mle.hdpm(hdpm.data$y, P=200)
theta_mle <- hdpm.aux.mle$theta_mle

hdpm.aux.mle.result <- hdpm.mle.result
hdpm.aux.mle.result['True','log-lik'] <- round(hdpm.aux.filter$loglik,3)
hdpm.aux.mle.result['MLE', 'log-lik'] <- round(hdpm.aux.mle$loglik,3)
hdpm.aux.mle.result['True',2:7] <- round(c(theta$D.phi0, theta$D.phi1, theta$I.phi1, theta$P.int, theta$D.var, theta$I.var),2)
hdpm.aux.mle.result['MLE', 2:7] <- round(c(theta_mle$D.phi0, theta_mle$D.phi1, theta_mle$I.phi1, theta_mle$P.int, theta_mle$D.var, theta_mle$I.var),2)
hdpm.aux.mle.result

if(save.results) {
    save(hdpm.particle.mle.result, file="../results/hdpm.particle.mle.result.Rda")
    save(hdpm.aux.mle.result,      file="../results/hdpm.aux.mle.result.Rda")
}
