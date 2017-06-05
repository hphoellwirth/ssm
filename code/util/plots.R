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
plot.weights <- function(weights, col='blue', xlab, ylab='weight') {
    
    plot(rowMeans(weights), type='l', col=col, xlab=xlab, ylab=ylab, xaxt="n", las=2)
    xticks <- axis(side=1, at=(1:nrow(weights)))
    abline(v=xticks , lty=3)
}


