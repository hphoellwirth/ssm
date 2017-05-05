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

# load state space models
source("models/localLevel.R")

# load filters
source("../filter/kalman.R")

# ----------------------------------------------------------------------
# Test Kalman filter on (univariate) local level model
# ----------------------------------------------------------------------

# generate local level data
llm.data <- gen.llm.data(n=100)

# use Kalman filter to estimate model states
llm.est <- kalman.filter(llm.data$y)

# plot result
plot(llm.data$y, type='l', col="red")
lines(llm.data$x, col="blue")
lines(llm.est$a, col="green")

# ----------------------------------------------------------------------
# Test Kalman filter on trivariate local level model
# ----------------------------------------------------------------------

# generate local level data
cov.eta <- matrix(c(1,0.5,0,0.5,1,0.1,0,0.1,2), nrow=3, ncol=3)
mllm.data <- gen.multi.llm.data(n=100, d=3, cov.eta=cov.eta)

# use Kalman filter to estimate model states
# llm.est <- kalman.filter(mllm.data$y)

# plot result
# plot(llm.data$y, type='l', col="red")
# lines(llm.data$x, col="blue")
# lines(llm.est$a, col="green")

