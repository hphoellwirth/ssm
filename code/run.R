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
source("models/01_localLevel.R")

# load filters
source("../filter/kalman.R")

# ----------------------------------------------------------------------
# Test Kalman filter on local level model
# ----------------------------------------------------------------------

# generate local level data
llm.data <- gen.llm.data(n=100)

# use Kalman filter to estimate model states
llm.est <- kalman.filter(llm.data$y)

# plot result
plot(llm.data$y, type='l', col="red")
lines(llm.data$alpha, col="blue")
lines(llm.est$a, col="green")

