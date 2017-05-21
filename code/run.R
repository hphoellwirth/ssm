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

# ----------------------------------------------------------------------
# Test Kalman filter on (univariate) local level model
# ----------------------------------------------------------------------

# generate local level data
T <- 100
var.eta <- 1.7
llm.data <- gen.llm.data(n=T, var.eta=var.eta)

# use Kalman filter to estimate model states
llm.kalman.filter <- kalman.filter(llm.data$y, cov.eta=var.eta)

# use particle filter to estimate model states
P <- 200
eta.sim <- matrix(rnorm(P*T, mean=0, sd=1), nrow=P, ncol=T) 
u.sim   <- matrix(runif(P*T, min=0, max=1), nrow=P, ncol=T)   
for (t in c(1:T)) {u.sim[,t] <- sort( u.sim[,t] )}
llm.particle.filter <- particle.filter(llm.data$y, cov.eta=var.eta, eta.sim=eta.sim, u.sim=u.sim)

# plot observations, states, and estimates
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(llm.data$y, type='l', col="red")
lines(llm.data$x, col="blue")
lines(llm.kalman.filter$a, col="green")
lines(llm.particle.filter$x.pr, col="orange")

# plot log-likelihood for different var.eta values
eta <- seq(0.5,2,0.1)
ll.kalman <- ll.particle <- rep(0,length(eta))
for (i in 1:length(eta)) {
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.particle[i] <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim)$loglik
    
}

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(eta, ll.kalman, type='b', col="green", xlab="var.eta with Kalman filter", ylab="log-likelihood", xaxt="n", las=2)
points(var.eta, llm.kalman.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=eta)
abline(v=xticks , lty=3)

plot(eta, ll.particle, type='b', col="orange", xlab="var.eta with particle filter", ylab="log-likelihood", xaxt="n", las=2)
points(var.eta, llm.particle.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=eta)
abline(v=xticks , lty=3)


# ----------------------------------------------------------------------
# Test Kalman filter on trivariate local level model
# ----------------------------------------------------------------------

# generate local level data
D <- 3
T <- 100
cov.eta.var <- c(4.2,2.8,0.9)
cov.eta.rho <- 0.7
mllm.data <- gen.multi.llm.data(n=T, d=D, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use Kalman filter to estimate model states
mllm.kalman.filter <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use particle filter to estimate model states
P <- 200
eta.sim <- u.sim <- list()
for (t in 1:T) {
    eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D)) 
    u.sim[[t]]   <- matrix(runif(P*D, min=0, max=1), nrow=P, ncol=D) 
    for (d in 1:D) {u.sim[[t]][,d] <- sort( u.sim[[t]][,d] )}
}
mllm.particle.filter <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim)
#mllm.particle.filter <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(1.0,2.8,4.5), 0.7), eta.sim=eta.sim, u.sim=u.sim)
#mllm.particle.filter <- m.particle.filter(mllm.data$y, cov.eta=diag(3), eta.sim=eta.sim, u.sim=u.sim)


# plot observations, states, and estimates
par(mfrow=c(3,1), mar=c(1,1,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red")
    lines(mllm.data$x[,d], col="blue")
    lines(mllm.kalman.filter$a[,d], col="green")
    lines(mllm.particle.filter$x.pr[,d], col="orange")
}

# plot log-likelihood for different rho values
rho <- seq(-1,1,0.1)
ll.kalman <- ll.particle <- rep(0,length(rho))
for (i in 1:length(rho)) {
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, rho[i]))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, rho[i]), eta.sim=eta.sim, u.sim=u.sim)$loglik
}

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(rho, ll.kalman, type='b', col="green", xlab="rho with Kalman filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.rho, mllm.kalman.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=rho)
abline(v=xticks , lty=3)
plot(rho, ll.particle, type='b', col="orange", xlab="rho with particle filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.rho, mllm.particle.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=rho)
abline(v=xticks , lty=3)

# plot log-likelihood for different eta.var1 values
var1 <- seq(2,5.6,0.2)
ll.kalman <- ll.particle <- rep(0,length(var1))
for (i in 1:length(var1)) {
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(var1[i],2.8,0.9), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(var1[i],2.8,0.9), cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim)$loglik
}

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(var1, ll.kalman, type='b', col="green", xlab="eta.var1 with Kalman filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.var[1], mllm.kalman.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=var1)
abline(v=xticks , lty=3)
plot(var1, ll.particle, type='b', col="orange", xlab="eta.var1 with particle filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.var[1], mllm.particle.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=var1)
abline(v=xticks , lty=3)

# plot log-likelihood for different eta.var3 values
var3 <- seq(0.2,4,0.2)
ll.kalman <- ll.particle <- rep(0,length(var3))
for (i in 1:length(var3)) {
    ll.kalman[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,2.8,var3[i]), cov.eta.rho))$loglik
    ll.particle[i] <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(c(4.2,2.8,var3[i]), cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim)$loglik
}

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(var3, ll.kalman, type='b', col="green", xlab="eta.var1 with Kalman filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.var[3], mllm.kalman.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=var3)
abline(v=xticks , lty=3)
plot(var3, ll.particle, type='b', col="orange", xlab="eta.var1 with particle filter", ylab="log-likelihood", xaxt="n", las=2)
points(cov.eta.var[3], mllm.particle.filter$loglik) # highlight true parameter
xticks <- axis(side=1, at=var3)
abline(v=xticks , lty=3)

# estimate model parameters, using Kalman filter
mllm.mle <- kalman.mle(mllm.data$y, 4)
print(paste('     True parameters:', round(cov.eta.var[1],2), round(cov.eta.var[2],2), round(cov.eta.var[3],2), round(cov.eta.rho,2)))
print(paste('Estimated parameters:', round(mllm.mle$theta_mle[1],2), round(mllm.mle$theta_mle[2],2), round(mllm.mle$theta_mle[3],2), round(mllm.mle$theta_mle[4],2)))
print(paste('     True log-likelihood:', round(mllm.kalman.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(mllm.mle$loglik,3)))

# estimate model parameters, using particle filter
mllm.mle <- particle.mle(mllm.data$y, 4, P)
print(paste('     True parameters:', round(cov.eta.var[1],2), round(cov.eta.var[2],2), round(cov.eta.var[3],2), round(cov.eta.rho,2)))
print(paste('Estimated parameters:', round(mllm.mle$theta_mle[1],2), round(mllm.mle$theta_mle[2],2), round(mllm.mle$theta_mle[3],2), round(mllm.mle$theta_mle[4],2)))
print(paste('     True log-likelihood:', round(mllm.kalman.filter$loglik,3)))
print(paste('Estimated log-likelihood:', round(mllm.mle$loglik,3)))

# plot filter with true and estimated parameters
mllm.est <- kalman.filter(mllm.data$y, cov.eta=construct.cov(mllm.mle$theta_mle[1:3], mllm.mle$theta_mle[4]))
par(mfrow=c(3,1), mar=c(1,1,1,1))
for (d in 1:3) {
    plot(mllm.data$x[,d], type='l', col="blue") 
    lines(mllm.filter$a[,d], col="green")
    lines(mllm.est$a[,d], col="purple")
}

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
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(as.vector(t(as.matrix(hdpm.data$x))), type='p', pch=19, col="red", xlab="time", ylab="counts / parameter (components)")
lines(as.vector(t(as.matrix(hdpm.data$lambda))), type='l', col="black")
lines(hdpm.data$D, type='l', col="orange")
lines(hdpm.data$P, type='l', col="green")
lines(hdpm.data$I, type='l', col="blue")

# plot log parameter (components)
plot(log(as.vector(t(as.matrix(hdpm.data$lambda)))), type='l', col="black", xlab="time", ylab="log parameter (components)")
lines(log(hdpm.data$D), type='l', col="orange")
lines(log(hdpm.data$P), type='l', col="green")
lines(log(hdpm.data$I), type='l', col="blue")


