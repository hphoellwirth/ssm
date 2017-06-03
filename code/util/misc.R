# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# Miscellaneous methods
#
# (Author) Hans-Peter Höllwirth
# (Date)   06.2017


# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

# load libraries
# NONE

# ----------------------------------------------------------------------
# Performance measure methods
# ----------------------------------------------------------------------
mse <- function(target, pred) {
    return(mean((target - pred)**2))
}

# ----------------------------------------------------------------------
# Matching methods
# ----------------------------------------------------------------------
match.rows <- function(tx, ty) {
    # matches rows in table tx with the (first occuring) row in table ty
    return(sapply(1:nrow(tx), function(r) which(apply(ty, 1, identical, tx[r,]))[1]))
}




