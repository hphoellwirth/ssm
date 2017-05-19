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
var.eta <- 1
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


# ----------------------------------------------------------------------
# Test Kalman filter on trivariate local level model
# ----------------------------------------------------------------------

# generate local level data
cov.eta.var <- c(4.3,2.8,0.9)
cov.eta.rho <- 0.3
mllm.data <- gen.multi.llm.data(n=100, d=3, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use Kalman filter to estimate model states
mllm.filter <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# plot observations, states, and estimates
par(mfrow=c(3,1), mar=c(1,1,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red")
    lines(mllm.data$x[,d], col="blue")
    lines(mllm.filter$a[,d], col="green")
}

# plot log-likelihood for different rho values
rho <- seq(-1,1,0.1)
ll <- rep(0,length(rho))
for (i in 1:length(rho)) {
    ll[i] <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, rho[i]))$loglik
}

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(rho, ll, type='b', col="red", xlab="rho", ylab="log-likelihood", xaxt="n", las=2)
xticks <- axis(side=1, at=rho)
abline(v=xticks , lty=3)

# estimate model parameters, using Kalman filter
mllm.mle <- mll.mle(mllm.data$y, 4)
print(paste('     True parameters:', round(cov.eta.var[1],2), round(cov.eta.var[2],2), round(cov.eta.var[3],2), round(cov.eta.rho,2)))
print(paste('Estimated parameters:', round(mllm.mle$theta_mle[1],2), round(mllm.mle$theta_mle[2],2), round(mllm.mle$theta_mle[3],2), round(mllm.mle$theta_mle[4],2)))
print(paste('     True log-likelihood:', round(mllm.filter$loglik,3)))
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


