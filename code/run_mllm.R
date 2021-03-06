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
save.results <- FALSE

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
mllm.data <- gen.multi.llm.data(T, d=D, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# plot observations and states (vertical)
if(save.plots) png("../images/mllm-realization.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(2,2,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red", ylim=c(min(mllm.data$y), max(mllm.data$y)))
    lines(mllm.data$x[,d], col="black")
    if (d == D) legend('bottomright', legend=c('observation','state'), col=c('red','black'), cex=1.0, lty=c(1,1), lwd=c(2,2))
}
if(save.plots) dev.off()

# plot observations and states (horizontal)
if(save.plots) png("../images/mllm-realization-v2.png", width=1000, height=300, pointsize=20)
par(mfrow=c(1,3), mar=c(2,2,1,1))
for (d in 1:D) {
    plot(mllm.data$y[,d], type='p', pch=19, col=alpha("red",0.7), main=paste('d',d), xlab="time", ylab=NULL, ylim=c(min(mllm.data$y), max(mllm.data$y)))
    lines(mllm.data$x[,d], col="black")
    if (d==1) legend('topleft', legend=c('observation','state'), col=c('red','black'), cex=1.0, lty=c(0,1), lwd=c(0,2.5), pch=c(19,NA))
}
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Use filters with true parameters to estimate model states
# ----------------------------------------------------------------------

# use Kalman filter to estimate model states
mllm.kalman.filter <- kalman.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho))

# use particle filter to estimate model states
P <- 50
eta.sim <- list()
for (t in 1:T) {eta.sim[[t]] <- mvrnorm(P, mu=rep(0,D), Sigma=diag(D))}
u.sim <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}
mllm.particle.filter <- m.particle.filter(mllm.data$y, cov.eta=construct.cov(cov.eta.var, cov.eta.rho), eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))

# plot observations, states, and predictions
if(save.plots) png("../images/mllm-predictions.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(2,2,1,1))
for (d in 1:3) {
    plot(mllm.data$y[,d], type='l', col="red", ylim=c(min(mllm.data$y), max(mllm.data$y)))
    lines(mllm.data$x[,d], col="black")
    lines(mllm.kalman.filter$x.pr[,d], col="green")
    lines(mllm.particle.filter$x.pr[,d], col="orange")
    if (d == D) legend('bottomright', legend=c('observation','state','Kalman filter', 'particle filter'), col=c('red','black','green','orange'), cex=1.0, lty=rep(1,4), lwd=rep(2,4))
}
if(save.plots) dev.off()

# plot Kalman predictions (with confidence interval)
if(save.plots) png("../images/mllm-estimate-kalman.png", width=1000, height=300, pointsize=20)
par(mfrow=c(1,3), mar=c(2,2,1,1))
for (d in 1:D) {
    plot(mllm.data$y[,d], type='p', pch=19, col=alpha("red",0.7), main=paste('d',d), xlab="time", ylab=NULL, ylim=c(min(mllm.data$y), max(mllm.data$y)))
    lines(mllm.particle.filter$x.up[,d], col='green', lwd=1.5)
    lines(mllm.kalman.filter$x.up[,d] + 2*sqrt(unlist(lapply(mllm.kalman.filter$x.up.cov, '[[', (1 + (d-1)*4)))), col='black', lty=2)
    lines(mllm.kalman.filter$x.up[,d] - 2*sqrt(unlist(lapply(mllm.kalman.filter$x.up.cov, '[[', (1 + (d-1)*4)))), col='black', lty=2)
    if (d==1) legend('topleft', legend=c('observation','estimate', '90% CI'), col=c('red','green','black'), cex=1.0, lty=c(0,1,2), lwd=c(0,2.5,2.5), pch=c(19,NA,NA))
}
if(save.plots) dev.off()

# plot SIR predictions (with confidence interval)
if(save.plots) png("../images/mllm-estimate-sir.png", width=1000, height=300, pointsize=20)
par(mfrow=c(1,3), mar=c(2,2,1,1))
for (d in 1:D) {
    states <- matrix(nrow=T, ncol=P)
    for (t in 1:(T)) {states[t,] <- mllm.particle.filter$x.up.particles[[t]][,d]}
    
    plot(mllm.data$y[,d], type='p', pch=19, col=alpha("red",0.7), main=paste('d',d), xlab="time", ylab=NULL, ylim=c(min(mllm.data$y), max(mllm.data$y)))
    lines(mllm.kalman.filter$x.up[,d], col='orange', lwd=1.5)
    lines(rowQuantiles(states, probs=c(0.95)), col='black', lty=2)
    lines(rowQuantiles(states, probs=c(0.05)), col='black', lty=2)
    if (d==1) legend('topleft', legend=c('observation','estimate', '90% CI'), col=c('red','orange','black'), cex=1.0, lty=c(0,1,2), lwd=c(0,2.5,2.5), pch=c(19,NA,NA))
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
set.seed(1000)
P <- 50
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

if(save.plots) png("../images/mllm-loglik-rho.png", width=500, height=400, pointsize=20)
par(mfrow=c(3,1), mar=c(2,2,1,1))
plot.loglik.zoom(rho, ll.kalman, cov.eta.rho, 'green', 'rho with Kalman filter', 'Kalman')
plot.loglik.zoom(rho, ll.particle, cov.eta.rho, 'orange', 'rho with particle filter','SIR')
plot.loglik.zoom(rho, ll.aux, cov.eta.rho, 'magenta', 'rho with auxiliary filter', 'IS particle')
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

if(save.plots) png("../images/mllm-loglik-var1.png", width=1000, height=750, pointsize=20)
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

if(save.plots) png("../images/mllm-loglik-var2.png", width=500, height=400, pointsize=20)
par(mfrow=c(3,1), mar=c(2,4,1,1))
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

if(save.plots) png("../images/mllm-loglik-var3.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(var3, ll.kalman, cov.eta.var[3], 'green', 'eta.var3 with Kalman filter')
plot.loglik(var3, ll.particle, cov.eta.var[3], 'orange', 'eta.var3 with particle filter')
plot.loglik(var3, ll.aux, cov.eta.var[3], 'magenta', 'eta.var3 with auxiliary filter')
if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Compute MLE with different filters
# ----------------------------------------------------------------------
mllm.mle.result <- data.frame(matrix(nrow=2, ncol=5), row.names=c('True','MLE'))
colnames(mllm.mle.result) <- c('log-lik','cov.eta.var1','cov.eta.var2','cov.eta.var3','cov.eta.rho')

# estimate model parameters, using Kalman filter
mllm.kalman.mle <- kalman.mle(mllm.data$y, D=D)

mllm.kalman.mle.result <- mllm.mle.result
mllm.kalman.mle.result['True','log-lik'] <- round(mllm.kalman.filter$loglik,3)
mllm.kalman.mle.result['MLE', 'log-lik'] <- round(mllm.kalman.mle$loglik,3)
mllm.kalman.mle.result['True',2:5] <- round(c(cov.eta.var[1], cov.eta.var[2], cov.eta.var[3], cov.eta.rho),2)
mllm.kalman.mle.result['MLE', 2:5] <- round(mllm.kalman.mle$theta_mle,2)
mllm.kalman.mle.result

# estimate model parameters, using particle filter
mllm.particle.mle <- m.particle.mle(mllm.data$y, D=D, P=50)

mllm.particle.mle.result <- mllm.mle.result
mllm.particle.mle.result['True','log-lik'] <- round(mllm.particle.filter$loglik,3)
mllm.particle.mle.result['MLE', 'log-lik'] <- round(mllm.particle.mle$loglik,3)
mllm.particle.mle.result['True',2:5] <- round(c(cov.eta.var[1], cov.eta.var[2], cov.eta.var[3], cov.eta.rho),2)
mllm.particle.mle.result['MLE', 2:5] <- round(mllm.particle.mle$theta_mle,2)
mllm.particle.mle.result

# estimate model parameters, using auxiliary filter
mllm.aux.mle <- m.aux.mle(mllm.data$y, D=D, P=50)

mllm.aux.mle.result <- mllm.mle.result
mllm.aux.mle.result['True','log-lik'] <- round(mllm.aux.filter$loglik,3)
mllm.aux.mle.result['MLE', 'log-lik'] <- round(mllm.aux.mle$loglik,3)
mllm.aux.mle.result['True',2:5] <- round(c(cov.eta.var[1], cov.eta.var[2], cov.eta.var[3], cov.eta.rho),2)
mllm.aux.mle.result['MLE', 2:5] <- round(mllm.aux.mle$theta_mle,2)
mllm.aux.mle.result

if(save.results) {
    save(mllm.kalman.mle.result,   file="../results/mllm.kalman.mle.result.Rda")
    save(mllm.particle.mle.result, file="../results/mllm.particle.mle.result.Rda")
    save(mllm.aux.mle.result,      file="../results/mllm.aux.mle.result.Rda")
}

# ----------------------------------------------------------------------
# Compare filters with true and estimated (MLE) parameters
# ----------------------------------------------------------------------

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





