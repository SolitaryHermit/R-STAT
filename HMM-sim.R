# Arguments of HMM Model
dice <- c("F", "L")
P <- matrix(c(0.99, 0.01, 0.25, 0.75), nrow=2, ncol=2, byrow=TRUE)
rownames(P) <- dice
colnames(P) <- dice
fair_prob <- rep(1/6, 6)
loaded_prob <- c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)
init_prob <- c(1/2, 1/2)

# 
markov <- function(init, P, n) {
    simstatus <- rep(0, n)
    simstatus[1] <- sample(0:4, 1, prob=init)
    for (i in 2:n) {
        simstatus[i] <- sample(0:4, 1, prob=P[simstatus[i-1]+1,])
    }
    return(simstatus)
}
