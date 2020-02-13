# Arguments of HMM Model
dice <- c("F", "L")
number <- seq(1, 6)
# Setting states and transition matrix
init_prob <- c(1/2, 1/2)
P <- matrix(c(0.99, 0.01, 0.25, 0.75), nrow=2, ncol=2, byrow=TRUE)
rownames(P) <- dice
colnames(P) <- dice
# Setting rolling probabilities
fair_roll <- rep(1/6, 6)
loaded_roll <- c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)
roll_prob <- matrix(c(fair_roll, loaded_roll), nrow=2, ncol=6, byrow=TRUE)
rownames(roll_prob) <- dice
colnames(roll_prob) <- number

# Implementation of HMM Model
hmm <- function(n, init_prob, trans_mat, outcome_mat, state_label, outcome_label) {
    # Set up results
    states <- rep(0, n)
    outcomes <- rep(0, n)
    # Initialize the states
    states[1] <- state_label[sample(1:length(state_label), 1, prob=init_prob)]
    outcomes[1] <- outcome_label[sample(1:length(outcome_label), 1, prob=outcome_mat[states[1]])]
    # Continue simulating
    for (i in 2:n) {
        states[i] <- state_label[sample(1:length(state_label), 1, prob=trans_mat[states[i-1]])]
        outcomes[i] <- outcome_label[sample(1:length(outcome_label), 1, prob=outcome_mat[states[i]])]
    }
    return(list(states, outcomes))
}
