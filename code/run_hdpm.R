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
source("filter/kalman.R")
source("filter/particle.R")
source("filter/aux.R")

# load utilities
source("util/plots.R")


# ----------------------------------------------------------------------
# Simulate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------

# generate hierarchical dynamic Poisson data
set.seed(1000)
hdpm.data <- gen.hdpm.data(N=5, M=20, D.phi0=0.7, D.phi1=0.6, I.phi1=0.3, P.int=0.8, D.var=1, I.var=1, a1=0, P1=1)

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



