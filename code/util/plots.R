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
plot.loglik <- function(para, loglik, true.para, true.loglik, col='blue', xlab='parameter') {
    # remove points further than 10% away from the median
    medianll <- median(loglik, na.rm=TRUE)   
    loglik <- replace(loglik, (loglik < 1.1*medianll | loglik > 0.9*medianll), NaN)
    
    plot(para, loglik, type='b', col=col, xlab=xlab, ylab="log-likelihood", xaxt="n", las=2)
    points(true.para, true.loglik) # highlight true parameter
    xticks <- axis(side=1, at=para)
    abline(v=xticks , lty=3)
}


