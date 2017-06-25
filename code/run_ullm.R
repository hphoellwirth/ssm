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
              upperCL=(llm.kalman.filter$x.up + 2*sqrt(unlist(llm.kalman.filter$x.up.cov))), 
              lowerCL=(llm.kalman.filter$x.up - 2*sqrt(unlist(llm.kalman.filter$x.up.cov))), 
              ylab='Kalman estimate', col='green')
if(save.plots) dev.off()

# use particle filter to estimate model states
P <- 200
llm.particle.filter <- particle.filter(llm.data$y, cov.eta=var.eta, P=P, x_up.init=rep(0,P), use.csir=FALSE)

if(save.plots) png("../images/ullm-estimate-particle.png", width=600, height=450, pointsize=14)
plot.estimate(llm.data$y, llm.particle.filter$x.up, 
              upperCL=rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.95)), 
              lowerCL=rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.05)),  
              ylab='SIR particle estimate', col='orange')
if(save.plots) dev.off()

# plot comparison of Kalman and particle filter
if(save.plots) png("../images/ullm-filter-comparison.png", width=600, height=450, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(llm.data$x, type='l', col="black", xlab="time", ylab='estimates')
lines(llm.kalman.filter$x.up, col='green')
lines(llm.particle.filter$x.up, col='orange')
lines(llm.kalman.filter$x.up + 2*sqrt(unlist(llm.kalman.filter$x.up.cov)), col=alpha('green',0.7), lty=2)
lines(llm.kalman.filter$x.up - 2*sqrt(unlist(llm.kalman.filter$x.up.cov)), col=alpha('green',0.7), lty=2)
lines(rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.95)), col=alpha('orange',0.7), lty=2)
lines(rowQuantiles(llm.particle.filter$x.up.particles, probs=c(0.05)), col=alpha('orange',0.7), lty=2)
legend('topleft', legend=c('latent states','Kalman estimate', 'Kalman 95% CL', 'SIR estimate', 'SIR 95% CL'), col=c('black','green','green','orange','orange'), cex=1.0, lty=c(1,1,2,1,2), lwd=2.5)
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
llm.aux.filter <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=var.eta, cov.eta.aux=1)

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



# ----------------------------------------------------------------------
# Produce log-likelihood plots to compare filters
# ----------------------------------------------------------------------

# draw standard normal transition noise and particles for comparison of particle filter
eta.sim <- matrix(rnorm(T*P, mean=0, sd=1), nrow=T, ncol=P) 
u.sim   <- matrix(runif(T*P, min=0, max=1), nrow=T, ncol=P)   
for (t in c(1:T)) {u.sim[t,] <- sort( u.sim[t,] )}

# compute log-likelihood for wide range of different var.eta values
eta <- seq(0.5,2,0.1)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.particle[i] <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P))$loglik
    ll.aux[i]      <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=var.eta)$loglik
}

# plot log-likelihood
if(save.plots) png("../images/univariate-loglik.png", width=1000, height=750, pointsize=20)
par(mfrow=c(3,1), mar=c(4,4,1,1))
plot.loglik(eta, ll.kalman, var.eta, 'green', 'var.eta with Kalman filter')
plot.loglik(eta, ll.particle, var.eta, 'orange', 'var.eta with particle filter')
plot.loglik(eta, ll.aux, var.eta, 'magenta', 'var.eta with auxiliary filter')
if(save.plots) dev.off()

# plot log-likelihood for different var.eta values zoomed around true value
eta <- seq(1.350,1.450,0.001)
ll.kalman <- ll.particle <- ll.aux <- rep(0,length(eta))
for (i in 1:length(eta)) {
    cat('.')
    ll.kalman[i]   <- kalman.filter(llm.data$y, cov.eta=eta[i])$loglik
    ll.particle[i] <- particle.filter(llm.data$y, cov.eta=eta[i], eta.sim=eta.sim, u.sim=u.sim, x_up.init=rep(0,P), use.csir = FALSE)$loglik 
    ll.aux[i]      <- aux.filter(llm.data$y, x.pr=llm.particle.filter$x.pr.particles, x.up=llm.particle.filter$x.up.particles, cov.eta=eta[i], cov.eta.aux=var.eta)$loglik
}

ll.kalman   <- ll.kalman / T
ll.particle <- ll.particle / T
ll.aux      <- ll.aux / T
ll.kalman   <- ll.kalman / abs( ll.kalman[ which(round(eta,3)==var.eta)] )
ll.particle <- ll.particle / abs( ll.particle[ which(round(eta,3)==var.eta)] )
ll.aux      <- ll.aux / abs( ll.aux[ which(round(eta,3)==var.eta)] )

# plot zoomed log-likelihood
if(save.plots) png("../images/univariate-loglik-detail.png", width=750, height=500, pointsize=15)
par(mfrow=c(1,1), mar=c(4,4,1,1))
matplot(eta, cbind(ll.kalman,ll.particle,ll.aux), type='b', col=c("green","orange","magenta") , ylab="log-likelihood", xaxt="n", las=2)
points(var.eta, -1 ) # highlight true parameter
xticks <- axis(side=1, at=eta)
abline(v=xticks , lty=3)
if(save.plots) dev.off()


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
Ts <- c(50,100,150,200,250,300,400,500)
Ps <- c(100,200,300,400,500)
R <- 10 # runs per combination

# generate R local level observations for each length T
var.eta <- 1.4
set.seed(1000)
llm.data <- list()
for (t in 1:length(Ts)) {
    llm.data[[t]] <- matrix(nrow=R, ncol=Ts[t])
    for (i in 1:R) {
        llm.data[[t]][i,] <- gen.llm.data(n=Ts[t], var.eta)$y
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


# compute MSE for particle filter (over parameters)
Ts <- c(50,100,150,200,250)
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





