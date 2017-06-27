# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Run space state methods on univariate local level model
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
library(matrixStats)
library(scales)

# if interactive, during the development, set to TRUE
interactive <- TRUE
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
# Generate univariate local level data
# ----------------------------------------------------------------------
set.seed(1000)
T <- 100
var.eta <- 1.4
llm.data <- gen.llm.data(T, var.eta=var.eta)

# plot system
if(save.plots) png("../images/ullm-realization.png", width=600, height=450, pointsize=14)
plot.realization(llm.data$y, llm.data$x)
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Use filters with true parameters to estimate model states
# ----------------------------------------------------------------------

# use Kalman filter to estimate model states
llm.kalman.filter <- kalman.filter(llm.data$y, cov.eta=var.eta)

if(save.plots) png("../images/ullm-estimate-kalman.png", width=600, height=450, pointsize=14)
plot.estimate(llm.data$y, llm.kalman.filter$x.up, 
              upper90CI=(llm.kalman.filter$x.up + 1.645*sqrt(unlist(llm.kalman.filter$x.up.cov))), 
              lower90CI=(llm.kalman.filter$x.up - 1.645*sqrt(unlist(llm.kalman.filter$x.up.cov))), 
              ylab='Kalman estimate', col='green')
if(save.plots) dev.off()

# use particle filter to estimate model states
P <- 200
llm.particle.filter <- particle.filter(llm.data$y, cov.eta=var.eta, P=P, x_up.init=rep(0,P), use.csir=FALSE)

if(save.plots) png("../images/ullm-estimate-particle.png", width=600, height=450, pointsize=14)
plot.estimate(llm.data$y, llm.particle.filter$x.up, 
              upper90CI=rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.95)), 
              lower90CI=rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.05)),  
              ylab='SIR particle estimate', col='orange')
if(save.plots) dev.off()

# plot comparison of Kalman and particle filter
if(save.plots) png("../images/ullm-filter-comparison.png", width=600, height=450, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(llm.data$x, type='l', col="black", xlab="time", ylab='estimates')
lines(llm.kalman.filter$x.up, col='green', lwd=1.5)
lines(llm.particle.filter$x.up, col='orange', lwd=1.5)
lines(llm.kalman.filter$x.up + 1.645*sqrt(unlist(llm.kalman.filter$x.up.cov)), col=alpha('green',0.7), lty=2)
lines(llm.kalman.filter$x.up - 1.645*sqrt(unlist(llm.kalman.filter$x.up.cov)), col=alpha('green',0.7), lty=2)
lines(rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.95)), col=alpha('orange',0.7), lty=2)
lines(rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.05)), col=alpha('orange',0.7), lty=2)
legend('topleft', legend=c('latent states','Kalman estimate', 'Kalman 90% CI', 'SIR estimate', 'SIR 90% CI'), col=c('black','green','green','orange','orange'), cex=1.0, lty=c(1,1,2,1,2), lwd=2.5)
if(save.plots) dev.off()

# plot observations, states, and predictions
if(save.plots) png("../images/ullm-predictions.png", width=1000, height=600, pointsize=14)
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(llm.data$y, type='p', pch=19, col="red")
lines(llm.data$x, col="black")
lines(llm.kalman.filter$x.pr, col="green")
lines(llm.particle.filter$x.pr, col="orange")
legend(10,15, c('observation','state','kalman','particle'), cex=1.0, lty=rep(1,4), lwd=rep(2.5,4), col=c('red','black','green','orange'))
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Evaluate the correctness of the auxiliary filter
# ----------------------------------------------------------------------

# use auxiliary filter to compute true log-likelihood
llm.aux.filter <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=var.eta, cov.eta.aux=1.0)

# plot importance weights over time
ylim <- c(min(rowMeans(llm.aux.filter$is.up), rowMeans(llm.aux.filter$is.pr)), max(rowMeans(llm.aux.filter$is.up), rowMeans(llm.aux.filter$is.pr)))
if(save.plots) png(paste0("../images/ullm_is_filt_weights_P",P,".png"), width=750, height=160, pointsize=14)
par(mfrow=c(1,1), mar=c(1,4,1,1))
plot.weights(llm.aux.filter$is.up, xlab='time t', ylab='filtering weights', ylim=ylim)
if(save.plots) dev.off()

if(save.plots) png(paste0("../images/ullm_is_pred_weights_P",P,".png"), width=750, height=200, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot.weights(llm.aux.filter$is.pr, xlab='time t', ylab='predictive weights', ylim=ylim)
if(save.plots) dev.off()

# observe change in log-likelihood as auxiliary parameter changes
eta <- seq(0.8,3,0.1)
ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.aux[i] <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=var.eta, cov.eta.aux=eta[i])$loglik
}

if(save.plots) png("../images/ullm-loglik-aux-eta.png", width=500, height=300, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot.loglik(eta, ll.aux, var.eta, 'magenta', 'var.eta.aux')
if(save.plots) dev.off()

# plot log-likelihood of auxiliary filter for different auxiliary parameters
eta.aux <- c(0.8,1.4,2.8)
eta <- seq(0.8,2,0.1)
ll.aux <- matrix(nrow=3, ncol=length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.aux[1,i]  <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=eta.aux[1])$loglik
    ll.aux[2,i]  <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=eta.aux[2])$loglik
    ll.aux[3,i]  <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=eta.aux[3])$loglik
}

if(save.plots) png("../images/ullm-loglik-aux-eta-2.png", width=500, height=300, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot.loglik.multi(eta, ll.aux, var.eta, col=c('red', 'magenta','brown'), 'var.eta', legend=eta.aux, legend.title='var.eta.aux')
if(save.plots) dev.off()

# ----------------------------------------------------------------------
# Produce log-likelihood plots to compare filters
# ----------------------------------------------------------------------

# draw standard normal transition noise and particles for comparison of particle filter
set.seed(1000)
eta.sim <- matrix(rnorm(T*P, mean=0, sd=1), nrow=T, ncol=P) 
u.sim   <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}

# compute log-likelihood for wide range of different var.eta values
eta <- seq(0.5,2,0.1)
ll.kalman <- ll.sir <- ll.csir <- ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.sir[i]      <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=FALSE)$loglik
    ll.csir[i]     <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=TRUE)$loglik
    ll.aux[i]      <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=var.eta)$loglik
}

# plot log-likelihood
if(save.plots) png("../images/ullm-loglik-eta.png", width=500, height=600, pointsize=20)
par(mfrow=c(4,1), mar=c(2,4,1,1))
plot.loglik(eta, ll.kalman, var.eta, 'green', 'var.eta with Kalman filter')
plot.loglik(eta, ll.sir, var.eta, 'orange', 'var.eta with SIR filter')
plot.loglik(eta, ll.csir, var.eta, 'blue', 'var.eta with CSIR filter')
plot.loglik(eta, ll.aux, var.eta, 'magenta', 'var.eta with IS particle filter')
if(save.plots) dev.off()

# plot log-likelihood for different var.eta values zoomed around true value
eta.det <- seq(1.350,1.450,0.001)
ll.kalman.zoom <- ll.sir.zoom <- ll.csir.zoom <- ll.aux.zoom <- rep(0,length(eta.det))
for (i in 1:length(eta.det)) {
    cat('.')
    ll.kalman.zoom[i] <- kalman.filter(llm.data$y, cov.eta=eta.det[i])$loglik
    ll.sir.zoom[i]    <- particle.filter(llm.data$y, cov.eta=eta.det[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=FALSE)$loglik 
    ll.csir.zoom[i]   <- particle.filter(llm.data$y, cov.eta=eta.det[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir=TRUE)$loglik 
    ll.aux.zoom[i]    <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta.det[i], cov.eta.aux=var.eta)$loglik
}

ll.kalman.zoom <- ll.kalman.zoom / T
ll.sir.zoom    <- ll.sir.zoom / T
ll.csir.zoom   <- ll.csir.zoom / T
ll.aux.zoom    <- ll.aux.zoom / T
ll.kalman.zoom <- ll.kalman.zoom / abs( ll.kalman.zoom[ which(round(eta.det,3)==var.eta)] )
ll.sir.zoom    <- ll.sir.zoom / abs( ll.sir.zoom[ which(round(eta.det,3)==var.eta)] )
ll.csir.zoom   <- ll.csir.zoom / abs( ll.csir.zoom[ which(round(eta.det,3)==var.eta)] )
ll.aux.zoom    <- ll.aux.zoom / abs( ll.aux.zoom[ which(round(eta.det,3)==var.eta)] )

# plot zoomed log-likelihood
if(save.plots) png("../images/ullm-loglik-zoom.png", width=500, height=600, pointsize=20)
par(mfrow=c(4,1), mar=c(2,1,1,1))
plot.loglik.zoom(eta.det, ll.kalman.zoom, var.eta, 'green', 'var.eta with Kalman filter', 'Kalman')
plot.loglik.zoom(eta.det, ll.sir.zoom, var.eta, 'orange', 'var.eta with SIR filter', 'SIR')
plot.loglik.zoom(eta.det, ll.csir.zoom, var.eta, 'blue', 'var.eta with CSIR filter', 'CSIR')
plot.loglik.zoom(eta.det, ll.aux.zoom, var.eta, 'magenta', 'var.eta with IS particle filter', 'IS particle')
if(save.plots) dev.off()

#if(save.plots) png("../images/ullm-loglik-detail.png", width=750, height=500, pointsize=15)
#par(mfrow=c(1,1), mar=c(4,4,1,1))
#matplot(eta, cbind(ll.kalman,ll.particle,ll.aux), type='b', col=c("green","orange","magenta") , ylab="log-likelihood", xaxt="n", las=2)
#points(var.eta, -1 ) # highlight true parameter
#xticks <- axis(side=1, at=eta)
#abline(v=xticks , lty=3)
#if(save.plots) dev.off()


# ----------------------------------------------------------------------
# Compute MLE with different filters
# ----------------------------------------------------------------------
llm.mle.result <- data.frame(matrix(nrow=2, ncol=2), row.names=c('True','MLE'))
colnames(llm.mle.result) <- c('log-lik','var.eta')

# estimate model parameter, using Kalman filter
llm.kalman.mle <- kalman.mle(llm.data$y, D=1)

llm.kalman.mle.result <- llm.mle.result
llm.kalman.mle.result['True','log-lik'] <- round(llm.kalman.filter$loglik,3)
llm.kalman.mle.result['MLE', 'log-lik'] <- round(llm.kalman.mle$loglik,3)
llm.kalman.mle.result['True','var.eta'] <- round(var.eta,2)
llm.kalman.mle.result['MLE', 'var.eta'] <- round(llm.kalman.mle$theta_mle,2)
llm.kalman.mle.result

# estimate model parameter, using particle filter
llm.particle.mle <- particle.mle(llm.data$y, P=100)

llm.particle.mle.result <- llm.mle.result
llm.particle.mle.result['True','log-lik'] <- round(llm.particle.filter$loglik,3)
llm.particle.mle.result['MLE', 'log-lik'] <- round(llm.particle.mle$loglik,3)
llm.particle.mle.result['True','var.eta'] <- round(var.eta,2)
llm.particle.mle.result['MLE', 'var.eta'] <- round(llm.particle.mle$theta_mle,2)
llm.particle.mle.result

# estimate model parameter, using auxiliary filter
llm.aux.mle <- aux.mle(llm.data$y, P=P)

llm.aux.mle.result <- llm.mle.result
llm.aux.mle.result['True','log-lik'] <- round(llm.aux.filter$loglik,3)
llm.aux.mle.result['MLE', 'log-lik'] <- round(llm.aux.mle$loglik,3)
llm.aux.mle.result['True','var.eta'] <- round(var.eta,2)
llm.aux.mle.result['MLE', 'var.eta'] <- round(llm.aux.mle$theta_mle,2)
llm.aux.mle.result

if(save.results) {
    save(llm.kalman.mle.result,   file="../results/ullm.kalman.mle.result.Rda")
    save(llm.particle.mle.result, file="../results/ullm.particle.mle.result.Rda")
    save(llm.aux.mle.result,      file="../results/ullm.aux.mle.result.Rda")
}

# ----------------------------------------------------------------------
# Compute filter MLE statistics for different T and P
# ----------------------------------------------------------------------
Ts <- c(50,100,250,500)
Ps <- c(20,50,200,500)
R <- 100 # runs per combination

# generate R local level observations for each length T
var.eta <- 1.4
set.seed(1000)
llm.data <- list()
for (t in 1:length(Ts)) {
    llm.data[[t]] <- matrix(nrow=R, ncol=Ts[t])
    for (i in 1:R) {
        llm.data[[t]][i,] <- gen.llm.data(T=Ts[t], var.eta)$y
    }
}

# compute MSE for Kalman filter (over parameters)
bias.kalman <- se.kalman <- mse.kalman <- rep(0,length(Ts))
for (t in 1:length(Ts)) {
    
    kalman.est <- rep(0,R)
    for (i in 1:R) {
        kalman.est[i] <- kalman.mle(llm.data[[t]][i,], D=1, verbose=FALSE)$theta_mle
    }
    
    bias.kalman[t] <- mean(kalman.est) - var.eta
    se.kalman[t] <- sd(kalman.est) / sqrt(R)
    mse.kalman[t] <- mse(kalman.est, rep(var.eta,R)) 
    cat('.')
} 

if(save.results) {
    save(bias.kalman, file="../results/ullm.kalman.bias.Rda")
    save(se.kalman,   file="../results/ullm.kalman.se.Rda")
    save(mse.kalman,  file="../results/ullm.kalman.mse.Rda")
}

if(save.plots) png("../images/ullm_mse_kalman.png", width=750, height=500, pointsize=15)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(Ts, mse.kalman, type='l', col='red', xlab='T', ylab='mse')
if(save.plots) dev.off()  


# compute MSE for CSIR particle filter (over parameters)
#Ts <- c(50,100,150,200,250)
bias.particle <- se.particle <- mse.particle <- data.frame(matrix(nrow=length(Ts), ncol=length(Ps)), row.names = paste0('T',Ts))
colnames(bias.particle) <- colnames(se.particle) <- colnames(mse.particle) <- paste0('P',Ps)
for (t in 1:length(Ts)) {
    for (p in 1:length(Ps)) {
        
        particle.est <- rep(0,R)
        for (i in 1:R) {
            particle.est[i] <- particle.mle(llm.data[[t]][i,], P=Ps[p], verbose=FALSE)$theta_mle
        }
        
        bias.particle[t,p] <- mean(particle.est) - var.eta
        se.particle[t,p] <- sd(particle.est) / sqrt(R)        
        mse.particle[t,p] <- mse(particle.est, rep(var.eta,R)) 
        cat('.')        
    }
} 

if(save.results) {
    save(bias.particle, file="../results/ullm.particle.bias.Rda")
    save(se.particle,   file="../results/ullm.particle.se.Rda")
    save(mse.particle,  file="../results/ullm.particle.mse.Rda")
}

# compute MSE for auxiliary filter (over parameters)
bias.aux <- se.aux <- mse.aux <- data.frame(matrix(nrow=length(Ts), ncol=length(Ps)), row.names = paste0('T',Ts))
colnames(bias.aux) <- colnames(se.aux) <- colnames(mse.aux) <- paste0('P',Ps)
for (t in 1:length(Ts)) {
    for (p in 1:length(Ps)) {
        
        aux.est <- rep(0,R)
        for (i in 1:R) {
            aux.est[i] <- aux.mle(llm.data[[t]][i,], P=Ps[p], verbose=FALSE)$theta_mle 
        }
        
        bias.aux[t,p] <- mean(aux.est) - var.eta
        se.aux[t,p] <- sd(aux.est) / sqrt(R)     
        mse.aux[t,p] <- mse(aux.est, rep(var.eta,R)) 
        cat('.')        
    }
} 

if(save.results) {
    save(bias.aux, file="../results/ullm.aux.bias.Rda")
    save(se.aux,   file="../results/ullm.aux.se.Rda")
    save(mse.aux,  file="../results/ullm.aux.mse.Rda")
}

# ----------------------------------------------------------------------
# Load and plot MLE statistics for different T and P
# ----------------------------------------------------------------------
load(file="../results/ullm.kalman.bias.Rda")
load(file="../results/ullm.kalman.se.Rda")
load(file="../results/ullm.kalman.mse.Rda")
load(file="../results/ullm.particle.bias.Rda")
load(file="../results/ullm.particle.se.Rda")
load(file="../results/ullm.particle.mse.Rda")
load(file="../results/ullm.aux.bias.Rda")
load(file="../results/ullm.aux.se.Rda")
load(file="../results/ullm.aux.mse.Rda")

if(save.plots) png("../images/ullm-mc-mle.png", width=1000, height=300, pointsize=20)
par(mfrow=c(1,3), mar=c(4,2,1,1))
# plot bias
bounds <- c(bias.kalman^2, bias.particle[,'P500']^2, bias.aux[,'P500']^2)
plot(Ts, bias.kalman^2, type='b', col='green', ylim=c(min(bounds),max(bounds)), xlab='T', ylab='', main='bias^2')
lines(Ts, bias.particle[,'P500']^2, type='b', col='blue', lwd='1.5')
lines(Ts, bias.particle[,'P200']^2, type='b', col=alpha('blue',0.6))
lines(Ts, bias.particle[,'P50']^2, type='b', col=alpha('blue',0.4))
lines(Ts, bias.particle[,'P20']^2, type='b', col=alpha('blue',0.2))
lines(Ts, bias.aux[,'P500']^2, type='b', col='magenta', lwd='1.5')
lines(Ts, bias.aux[,'P200']^2, type='b', col=alpha('magenta',0.6))
lines(Ts, bias.aux[,'P50']^2, type='b', col=alpha('magenta',0.4))
lines(Ts, bias.aux[,'P20']^2, type='b', col=alpha('magenta',0.2))
lines(Ts, bias.kalman^2, type='b', col='green', lwd='1.5')

# plot standard error
bounds <- c(se.kalman, se.particle[,'P500'], se.aux[,'P500'])
plot(Ts, se.kalman, type='b', col='green', ylim=c(min(bounds),max(bounds)), xlab='T', ylab='', main='standard error')
lines(Ts, se.particle[,'P500'], type='b', col='blue', lwd='1.5')
lines(Ts, se.particle[,'P200'], type='b', col=alpha('blue',0.6))
lines(Ts, se.particle[,'P50'], type='b', col=alpha('blue',0.4))
lines(Ts, se.particle[,'P20'], type='b', col=alpha('blue',0.2))
lines(Ts, se.aux[,'P500'], type='b', col='magenta', lwd='1.5')
lines(Ts, se.aux[,'P200'], type='b', col=alpha('magenta',0.6))
lines(Ts, se.aux[,'P50'], type='b', col=alpha('magenta',0.4))
lines(Ts, se.aux[,'P20'], type='b', col=alpha('magenta',0.2))
lines(Ts, se.kalman, type='b', col='green', lwd='1.5')

# plot MSE
bounds <- c(mse.kalman, mse.particle[,'P500'], mse.aux[,'P500'])
plot(Ts, mse.kalman, type='b', col='green', ylim=c(min(bounds),max(bounds)), xlab='T', ylab='', main='MSE')
lines(Ts, mse.particle[,'P500'], type='b', col='blue', lwd='1.5')
lines(Ts, mse.particle[,'P200'], type='b', col=alpha('blue',0.6))
lines(Ts, mse.particle[,'P50'], type='b', col=alpha('blue',0.4))
lines(Ts, mse.particle[,'P20'], type='b', col=alpha('blue',0.2))
lines(Ts, mse.aux[,'P500'], type='b', col='magenta', lwd='1.5')
lines(Ts, mse.aux[,'P200'], type='b', col=alpha('magenta',0.6))
lines(Ts, mse.aux[,'P50'], type='b', col=alpha('magenta',0.4))
lines(Ts, mse.aux[,'P20'], type='b', col=alpha('magenta',0.2))
lines(Ts, mse.kalman, type='b', col='green', lwd='1.5')
legend('topright',c('Kalman','CSIR (P=500)','IS (P=500)'), col=c('green','blue','magenta'), cex=1.0, lty=1, lwd=2.5)
if(save.plots) dev.off()



