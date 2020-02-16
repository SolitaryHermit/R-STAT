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

#------------------------------------------------------------------------------
# Description
#   Simulation of Hidden Markov Model
# Arguments:
#   n: Length of observations
#   init_prob: The initial probability
#   trans_mat: The transition matrix P
#   outcome_mat: The observation matrix E
#------------------------------------------------------------------------------
hidden_markov <- function(n, init_prob, trans_mat, outcome_mat) {
    # Store some quantities for convenience
    state_label <- colnames(trans_mat)   # The label of states, e.g., {"H","T"}
    outcome_label <- colnames(outcome_mat)   # The label of outcomes, e.g., {1,2,3,4,5,6}
    num_state <- length(state_label)    # The number of states
    num_outcome <- length(outcome_label)    # The number of outcomes
    # Set up objects to store the result
    states <- rep(0, n)
    obs <- rep(0, n)
    # Initialize the simulation
    states[1] <- state_label[sample(1:num_state, 1, prob=init_prob)]
    obs[1] <- outcome_label[sample(1:num_outcome, 1, prob=outcome_mat[states[1],])]
    # Continue simulating
    for (i in 2:n) {
        states[i] <- state_label[sample(1:num_state, 1, prob=trans_mat[states[i-1],])]
        obs[i] <- outcome_label[sample(1:num_outcome, 1, prob=outcome_mat[states[i],])]
    }
    return(list(Q=states, O=obs))
}

#------------------------------------------------------------------------------
# Description
#   Perform forward algorithm (alpha-pass) to a given sequence of 
#   observations to compute its probability
# Arguments:
#   obs: Sequence of observations
#   init_prob: The initial probability
#   trans_mat: The transition matrix P
#   outcome_mat: The observation matrix E
#------------------------------------------------------------------------------
forward_sum <- function(obs, init_prob, trans_mat, outcome_mat) {
    # Store some quantities for convenience
    len <- length(obs)    # Length of observations
    state_label <- colnames(trans_mat)   # The label of states, e.g., {"H","T"}
    num_state <- length(state_label)    # The number of states
    # Set up objects to store the results
    probs <- matrix(0, nrow=len, ncol=num_state)
    colnames(probs) <- state_label
    # Calculation for the first step
    for (k in 1:num_state) {
        probs[1,k] <- init_prob[k] * outcome_mat[k,obs[1]]
    }
    # Forward Summation
    for (i in 2:len) {
        for (current_state in state_label) {
            for (past_state in state_label) {
                probs[i,current_state] <- probs[i,current_state] + trans_mat[past_state,current_state] * probs[i-1,past_state]
            }
            probs[i,current_state] <- probs[i,current_state] * outcome_mat[current_state,obs[i]]
        }
    }
}
