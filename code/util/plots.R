# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Plotting methods
#
# (Author) Hans-Peter Höllwirth
# (Date)   22.05.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
# NONE

# ----------------------------------------------------------------------
# State/observation plots
# ----------------------------------------------------------------------
plot.realization <- function(observations, states, ylab='values') {
    par(mfrow=c(1,1), mar=c(4,4,1,1))
    plot(observations, type='p', pch=19, col="red", xlab="time", ylab=ylab)
    lines(states, col="black")
    legend('topleft', legend=c('observation','state'), col=c('red','black'), cex=1.0, lty=c(0,1), lwd=c(0,2.5), pch=c(19,NA))
}

plot.estimate <- function(observations, estimate, lowerCL, upperCL, ylab='estimate', col='blue') {
    par(mfrow=c(1,1), mar=c(4,4,1,1))
    plot(observations, type='p', pch=19, col="red", xlab="time", ylab=ylab)
    lines(estimate, col=col)
    lines(upperCL, col='black', lty=2)
    lines(lowerCL, col='black', lty=2)
    legend('topleft', legend=c('observation','estimate', 'bounds'), col=c('red',col,'black'), cex=1.0, lty=c(0,1,2), lwd=c(0,2.5,2.5), pch=c(19,NA,NA))
}

# ----------------------------------------------------------------------
# Likelihood plots
# ----------------------------------------------------------------------
plot.loglik <- function(para, loglik, true.para, col='blue', xlab='parameter') {
    # remove points further than 10% below the median
    medianll <- median(loglik, na.rm=TRUE)   
    loglik <- replace(loglik, (loglik < 1.1*medianll), NaN)
    
    plot(para, loglik, type='b', col=col, xlab=xlab, ylab="log-likelihood", xaxt="n", las=2)
    elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
    points(true.para, loglik[elementwise.all.equal(para, true.para)]) # highlight true parameter
    xticks <- axis(side=1, at=para)
    abline(v=xticks , lty=3)
}

# ----------------------------------------------------------------------
# Importance weight plot
# ----------------------------------------------------------------------
plot.weights <- function(weights, col='blue', xlab, ylab='weight', ylim=NA) {
    
    plot(rowMeans(weights), type='l', col=col, xlab=xlab, ylab=ylab, ylim=ylim, xaxt="n", las=2)
    xticks <- axis(side=1, at=(1:nrow(weights)))
    abline(v=xticks , lty=3)
}


