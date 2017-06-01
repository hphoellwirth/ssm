# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Run space state methods
#
# (Author) Hans-Peter Höllwirth
# (Date)   13.04.2017


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

# load models
source("models/localLevel.R")
source("models/hierarchDynPoisson.R")

# load filters
source("filter/kalman.R")
source("filter/particle.R")
source("filter/aux.R")

# load utilities
source("util/plots.R")

# ----------------------------------------------------------------------
# Test Kalman filter on (univariate) local level model
# ----------------------------------------------------------------------

# generate local level data
set.seed(1000)
T <- 100
var.eta <- 1.4
llm.data <- gen.llm.data(n=T, var.eta=var.eta)

# use Kalman filter to estimate model states
llm.kalman.filter <- kalman.filter(llm.data$y, cov.eta=var.eta)

# use particle filter to estimate model states
P <- 200
eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=P, ncol=T) 
u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
llm.particle.filter <- particle.filter(llm.data$y, cov.eta=var.eta, eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))

# use auxiliary filter to compute true log-likelihood
llm.aux.filter <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=var.eta, cov.eta.aux=var.eta)

# plot observations, states, and estimates
if(save.plots) png("../images/univariate-local-level.png", width=1000, height=600, pointsize=14)
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(llm.data$y, type='l', col="red")
lines(llm.data$x, col="blue")
lines(llm.kalman.filter$a, col="green")
lines(llm.particle.filter$x.pr, col="orange")
#lines(llm.aux.filter$x.pr, col="magenta")
legend(10,15, c('observation','state','kalman','particle','auxilary'), cex=1.0, lty=rep(1,4), lwd=rep(2.5,4), col=c('red','blue','green','orange','magenta'))
if(save.plots) dev.off()

# plot log-likelihood for different var.eta values
eta <- seq(0.5,2,0.1)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.particle[i] <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=var.eta)$loglik
}

if(save.plots) png("../images/univariate-loglik.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(eta, ll.kalman, var.eta, llm.kalman.filter$loglik, 'green', 'var.eta with Kalman filter')
plot.loglik(eta, ll.particle, var.eta, llm.particle.filter$loglik, 'orange', 'var.eta with particle filter')
plot.loglik(eta, ll.aux, var.eta, llm.aux.filter$loglik, 'magenta', 'var.eta with auxiliary filter')
if(save.plots) dev.off()

# plot log-likelihood for different var.eta values zoomed around true value
eta <- seq(1.350,1.450,0.001)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.particle[i] <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik 
    ll.aux[i]      <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=var.eta)$loglik
    
}

ll.kalman   <- ll.kalman / T
ll.particle <- ll.particle / T
ll.aux      <- ll.aux / T
ll.kalman   <- ll.kalman / abs( ll.kalman[ which(round(eta,3)==var.eta)] )
ll.particle <- ll.particle / abs( ll.particle[ which(round(eta,3)==var.eta)] )
ll.aux      <- ll.aux / abs( ll.aux[ which(round(eta,3)==var.eta)] )

if(save.plots) png("../images/univariate-loglik-detail.png", width=750, height=500, pointsize=15)
par(mfrow=c(1,1), mar=c(4,4,1,1))
matplot(eta, cbind(ll.kalman,ll.particle,ll.aux), type='b', col=c("green","orange","magenta") , ylab="log-likelihood", xaxt="n", las=2)
points(var.eta, -1 ) # highlight true parameter
xticks <- axis(side=1, at=eta)
abline(v=xticks , lty=3)
if(save.plots) dev.off()

# estimate model parameter, using Kalman filter
llm.kalman.mle <- kalman.mle(llm.data$y, 1)
print(paste('     True parameters:', var.eta))
print(paste('Estimated parameters:', round(llm.kalman.mle$theta_mle[1],3)))
print(paste('     True log-likelihood:', round(llm.kalman.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(llm.kalman.mle$loglik,3)))

# estimate model parameter, using particle filter
llm.particle.mle <- particle.mle(llm.data$y, P=P)
print(paste('     True parameters:', var.eta))
print(paste('Estimated parameters:', round(llm.particle.mle$theta_mle[1],3)))
print(paste('     True log-likelihood:', round(llm.particle.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(llm.particle.mle$loglik,3)))

# estimate model parameter, using auxiliary filter
llm.aux.mle <- aux.mle(llm.data$y, P=P)
print(paste('     True parameters:', var.eta))
print(paste('Estimated parameters:', round(llm.aux.mle$theta_mle[1],3)))
print(paste('     True log-likelihood:', round(llm.aux.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(llm.aux.mle$loglik,3)))

# TBD: compare MLE standard errors and execution time

# ----------------------------------------------------------------------
# Test Kalman filter on trivariate local level model
# ----------------------------------------------------------------------

# generate local level data
set.seed(1000)
D <- 3
T <- 100
cov.eta.var <- c(4.2,2.8,0.9)
cov.eta.rho <- 0.7
mllm.data <- gen.multi.llm.data(n=T, d=D, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use Kalman filter to estimate model states
mllm.kalman.filter <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use particle filter to estimate model states
P <- 200
eta.sim <- list()
for (t in 1:T) {eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D))}
u.sim <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
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

# plot log-likelihood for different rho values
rho <- seq(-1.0,1,0.1)
ll.kalman <- ll.particle <- rep(0,length(rho))
for (i in 1:length(rho)) {
    cat('.')
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, rho[i]))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, rho[i]), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-rho.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(rho, ll.kalman, cov.eta.rho, mllm.kalman.filter$loglik, 'green', 'rho with Kalman filter')
plot.loglik(rho, ll.particle, cov.eta.rho, mllm.particle.filter$loglik, 'orange', 'rho with particle filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var1 values
var1 <- seq(2,5.6,0.2)
ll.kalman <- ll.particle <- rep(0,length(var1))
for (i in 1:length(var1)) {
    cat('.')
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(var1[i],2.8,0.9), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(var1[i],2.8,0.9), cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var1.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(var1, ll.kalman, cov.eta.var[1], mllm.kalman.filter$loglik, 'green', 'eta.var1 with Kalman filter')
plot.loglik(var1, ll.particle, cov.eta.var[1], mllm.particle.filter$loglik, 'orange', 'eta.var1 with particle filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var2 values
var2 <- seq(1.6,4.6,0.2)
ll.kalman <- ll.particle <- rep(0,length(var2))
for (i in 1:length(var2)) {
    cat('.')
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,var2[i],0.9), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,var2[i],0.9), cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var2.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(var2, ll.kalman, cov.eta.var[2], mllm.kalman.filter$loglik, 'green', 'eta.var2 with Kalman filter')
plot.loglik(var2, ll.particle, cov.eta.var[2], mllm.particle.filter$loglik, 'orange', 'eta.var2 with particle filter')
if(save.plots) dev.off()

# plot log-likelihood for different eta.var3 values
var3 <- seq(0.2,4,0.2)
ll.kalman <- ll.particle <- rep(0,length(var3))
for (i in 1:length(var3)) {
    cat('.')
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,2.8,var3[i]), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,2.8,var3[i]), cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
}

if(save.plots) png("../images/trivariate-local-level-loglik-var3.png", width=1000, height=500, pointsize=14)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot.loglik(var3, ll.kalman, cov.eta.var[3], mllm.kalman.filter$loglik, 'green', 'eta.var3 with Kalman filter')
plot.loglik(var3, ll.particle, cov.eta.var[3], mllm.particle.filter$loglik, 'orange', 'eta.var3 with particle filter')
if(save.plots) dev.off()

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


# ----------------------------------------------------------------------
# Simulate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------

# generate hierarchical dynamic Poisson data
set.seed(1000)
hdpm.data <- gen.hdpm.data(n=5, m=20, D.phi0=0.7, D.phi1=0.6, I.phi1=0.3, P.int=0.8, D.var=1, I.var=1, a1=0, P1=1)

# plot data and parameter (components)
if(save.plots) png("../images/dyn-poisson.png", width=1000, height=500, pointsize=14)
par(mfrow=c(1,1), mar=c(2,2,1,1))
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



