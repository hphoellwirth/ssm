# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Hierarchical dynamic Poisson model
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017


# ----------------------------------------------------------------------
# Generate hierarchical dynamic Poisson model
# ----------------------------------------------------------------------
gen.hdpm.data <- function(N, M, theta, a1=0, P1=1) {
    log.lambda.D <- rep(0,N)
    log.lambda.P <- rep(0,M)
    log.lambda.I <- rep(0,N*M)
    
    # draw disturbances
    u.D <- rnorm(N, mean=0, sd=sqrt(theta$D.var))
    u.I <- rnorm(N*M, mean=0, sd=sqrt(theta$I.var))

    # draw intial state
    log.lambda.D[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    log.lambda.I[1] <- rnorm(1, mean=a1, sd=sqrt(P1))    
        
    # compute daily component
    for (n in 1:(N-1)) {
        log.lambda.D[n+1] <- theta$D.phi0 + theta$D.phi1 * log.lambda.D[n] + u.D[n+1] 
    }
    
    # compute periodic component
    for (m in 1:M) {
        log.lambda.P[m] <- sin(pi*(m-1)/M) * theta$P.int
    } 
    
    # compute intra-daily component
    for (nm in 1:(N*M-1)) {
        log.lambda.I[nm+1] <- theta$I.phi1 * log.lambda.I[nm] + u.I[nm+1] 
    }
    
    # compute states
    x <- data.frame(matrix(ncol = M, nrow = N))   
    log.lambda <- data.frame(matrix(ncol = M, nrow = N)) 
    for (n in 1:N) {
        for (m in 1:M) {
            log.lambda[n,m] <- log.lambda.D[n] + log.lambda.P[m] + log.lambda.I[(n-1)*M+m] # state
            x[n,m] <- rpois(1, lambda=exp(log.lambda[n,m])) # observation
        }
    }
    return(list(x=x, lambda=exp(log.lambda), D=exp(rep(log.lambda.D,each=M)), P=exp(rep(log.lambda.P,N)), I=exp(log.lambda.I)))
}




