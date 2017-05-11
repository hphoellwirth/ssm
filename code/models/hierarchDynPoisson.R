# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Hierarchical dynamic Poisson model
#
# (Author) Hans-Peter Höllwirth
# (Date)   28.04.2017


# ----------------------------------------------------------------------
# Generate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------
gen.hdpm.data <- function(n, m, D.phi0=0, D.phi1=0.5, I.phi1=0.5, P.int=1, D.var=1, I.var=1, a1=0, P1=1) {
    log.lambda.D <- rep(0,n)
    log.lambda.P <- rep(0,m)
    log.lambda.I <- rep(0,n*m)
    
    # draw intial state
    log.lambda.D[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    log.lambda.I[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    u.D <- rnorm(n, mean=0, sd=sqrt(D.var))
    u.I <- rnorm(n*m, mean=0, sd=sqrt(I.var))
    
    # compute daily component
    for (t in 1:(n-1)) {
        log.lambda.D[t+1] <- D.phi0 + D.phi1 * log.lambda.D[t] + u.D[t+1] 
    }
    
    # compute periodic component
    for (i in 1:m) {
        log.lambda.P[i] <- sin(pi*(i-1)/m) * P.int
    } 
    
    # compute intra-daily component
    for (ti in 1:(n*m-1)) {
        log.lambda.I[ti+1] <- I.phi1 * log.lambda.I[ti] + u.I[ti+1] 
    }
    
    # compute states
    x <- data.frame(matrix(ncol = m, nrow = n))   
    log.lambda <- data.frame(matrix(ncol = m, nrow = n)) 
    for (t in 1:n) {
        for (i in 1:m) {
            log.lambda[t,i] <- log.lambda.D[t] + log.lambda.P[i] + log.lambda.I[(t-1)*m+i] # parameter
            x[t,i] <- rpois(1, lambda=exp(log.lambda[t,i])) # state
        }
    }
    return(list(x=x, lambda=exp(log.lambda), D=exp(rep(log.lambda.D,each=m)), P=exp(rep(log.lambda.P,n)), I=exp(log.lambda.I)))
}




