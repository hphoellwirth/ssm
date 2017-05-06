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
interactive <- TRUE
if (interactive) {
    setwd("/Users/Hans-Peter/Documents/Masters/14D000/code")
} 

# load models
source("models/localLevel.R")
source("hierarchDynPoisson.R")

# load filters
source("../filter/kalman.R")

# ----------------------------------------------------------------------
# Test Kalman filter on (univariate) local level model
# ----------------------------------------------------------------------

# generate local level data
llm.data <- gen.llm.data(n=100)

# use Kalman filter to estimate model states
llm.est <- m.kalman.filter(llm.data$y)

# plot result
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(llm.data$y, type='l', col="red")
lines(llm.data$x, col="blue")
lines(llm.est$a, col="green")


# ----------------------------------------------------------------------
# Test Kalman filter on trivariate local level model
# ----------------------------------------------------------------------

# generate local level data
cov.eta <- matrix(c(1,0.95,0,0.95,1,0.1,0,0.1,2), nrow=3, ncol=3)
mllm.data <- gen.multi.llm.data(n=100, d=3, cov.eta=cov.eta)

# use Kalman filter to estimate model states
mllm.est <- m.kalman.filter(mllm.data$y)

# plot result
par(mfrow=c(3,1), mar=c(1,1,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red")
    lines(mllm.data$x[,d], col="blue")
    lines(mllm.est$a[,d], col="green")
}

# ----------------------------------------------------------------------
# Simulate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------

# generate hierarchical dynamic Poisson data
set.seed(1000)
hdpm.data <- gen.hdpm.data(n=5, m=20, D.phi0=0.7, D.phi1=0.6, I.phi1=0.3, P.int=0.8, D.var=1, I.var=1, a1=0, P1=1)

#x <- hdpm.data$x
#lambda <- hdpm.data$lambda
#D <- hdpm.data$D
#P <- hdpm.data$P
#I <- hdpm.data$I

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


