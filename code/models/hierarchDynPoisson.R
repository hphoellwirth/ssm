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
gen.hdpm.data <- function(N, M, D.phi0=0, D.phi1=0.5, I.phi1=0.5, P.int=1, D.var=1, I.var=1, a1=0, P1=1) {
    log.lambda.D <- rep(0,N)
    log.lambda.P <- rep(0,M)
    log.lambda.I <- rep(0,N*M)
    
    # draw intial state
    log.lambda.D[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    log.lambda.I[1] <- rnorm(1, mean=a1, sd=sqrt(P1))
    
    # draw disturbances
    u.D <- rnorm(N, mean=0, sd=sqrt(D.var))
    u.I <- rnorm(N*M, mean=0, sd=sqrt(I.var))
    
    # compute daily component
    for (n in 1:(N-1)) {
        log.lambda.D[n+1] <- D.phi0 + D.phi1 * log.lambda.D[n] + u.D[n+1] 
    }
    
    # compute periodic component
    for (m in 1:M) {
        log.lambda.P[m] <- sin(pi*(m-1)/M) * P.int
    } 
    
    # compute intra-daily component
    for (nm in 1:(N*M-1)) {
        log.lambda.I[nm+1] <- I.phi1 * log.lambda.I[nm] + u.I[nm+1] 
    }
    
    # compute states
    x <- data.frame(matrix(ncol = M, nrow = N))   
    log.lambda <- data.frame(matrix(ncol = M, nrow = N)) 
    for (n in 1:N) {
        for (m in 1:M) {
            log.lambda[n,m] <- log.lambda.D[n] + log.lambda.P[m] + log.lambda.I[(n-1)*M+m] # parameter
            x[n,m] <- rpois(1, lambda=exp(log.lambda[n,m])) # state
        }
    }
    return(list(x=x, lambda=exp(log.lambda), D=exp(rep(log.lambda.D,each=M)), P=exp(rep(log.lambda.P,N)), I=exp(log.lambda.I)))
}




