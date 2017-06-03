# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Sequential importance resampling algorithms
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------
dyn.load("sampler/csir.so")

# ----------------------------------------------------------------------
# Univariate sequential importance resampling algorithm
# ----------------------------------------------------------------------
sir <- function(p, w, u) {
    # filter particles from set of particles (p), 
    # using their importance weights (w) and uniformly drawn samples (u)
    P <- length(p)
    s <- rep(0,P)  
    
    # sorting and weighting (leads to smoother log-likelihood curve)
    w <- w / sum(w) # normalize weights
    p.sort <- cbind(seq(1,P,1),c(p)) 
    p.idx  <- p.sort[order(p.sort[,2]),]
    p <- p.idx[,2]
    w <- w[p.idx[,1]]
    w.cum <- c(0, cumsum(w))
    
    # importance resampling of particles
    j <- 1
    for (i in 1:P){
        while((w.cum[i] < u[j]) && (u[j] <= w.cum[i+1])){
            s[j] <- p[i]
            if (j < P) {
                j <- j + 1
            }
            else break
        }
    }
    return(s) # return filtered particles
}

# ----------------------------------------------------------------------
# Multivariate sequential importance resampling algorithm
# ----------------------------------------------------------------------
m.sir <- function(p, w, u) {
    # filter (multivariate) particles from set of particles (p), 
    # using their importance weights (w) and uniformly drawn samples (u)
    D <- ncol(p)
    P <- nrow(p)
    s <- matrix(0,P,D)
    
    # compute cumulative normalized weights
    w <- w / sum(w) 
    w.cum <- c(0, cumsum(w))
    
    # importance resampling of particles
    j <- 1
    for (i in 1:P){
        while((w.cum[i] < u[j]) && (u[j] <= w.cum[i+1])){
            s[j,] <- p[i,]
            if (j < P) {
                j <- j + 1
            }
            else break
        }
    }
    return(s) # return filtered particles
}

# ----------------------------------------------------------------------
# Continous (univariate) sequential importance resampling algorithm
# ----------------------------------------------------------------------
csir <- function(p, w, u) {
    P <- length(p)
    s <- rep(0,P)
    
    # sorting and weighting
    w <- w / sum(w)
    p.sort <- cbind(seq(1,P,1),p)
    p.idx <- p.sort[order(p.sort[,2]),]
    p <- p.idx[,2]
    alpha_idx <- 
    w <- w[p.idx[,1]]
    w.cum <- c(0, cumsum(w))
    p  <- c(p[1], p)
    
    # importance resampling of particles
    j <- 1
    for (i in 1:P){
        while((w.cum[i] < u[j]) & (u[j] <= w.cum[i+1])){
            s[j] <- p[i] + ((p[i+1]-p[i])/(w.cum[i+1]-w.cum[i])) * (u[j]-w.cum[i])
            if (j < P) {
                j <- j + 1
            }
            else break
        }
    }
    return(s)
}

# C implementation (crashes for low P<120 values)
csir.c <- function(p, w, u) {
    P <- length(p)
    s <- rep(0,P)
    
    # call C function
    results <- .C("csir", s=as.double(s), p=as.double(p), w=as.double(w), u=as.double(u), len=as.integer(P))
    return (results$s)     
}

