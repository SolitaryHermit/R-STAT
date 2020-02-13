#------------------------------------------------------------
# NAME: bindecomp
# UTILITY: Transform a decimal number into a binary vector
# ARGUMENT: n -- decimal number
#------------------------------------------------------------
bindecomp <- function(n) {
    if (n <= 0) {
        warning("The number should be positive!")
        return()
    }
    bin <- rep(0, as.integer(log2(n)) + 1)
    while (n > 0) {
        index <- as.integer(log2(n))
        bin[index + 1] <- 1
        n <- n - 2^index
    }
    return(bin)
}

#------------------------------------------------------------
# NAME: matpow
# UTILITY: Calculate the n-th power of some square matrix
# ARGUMENT: M -- matrix; n -- exponent
#------------------------------------------------------------
matpow <- function(M, n) {
    # Check whether M is a square matrix 
    if (nrow(M) != ncol(M)) {
        warning("The size of input matrix is invalid!")
        return()
    } 
    bin <- bindecomp(n) # Obtain the binary expansion of n
    power <- diag(1, nrow(M), ncol(M))  # Initialize power as identity matrix
    for (i in 1:length(bin)) {
        if (bin[i] == 1) {
          power <- power %*% M
        }
        M = M %*% M
    }
    return(power)
}
