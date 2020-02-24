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
#   outcome_mat: The emission matrix E
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
#   Perform forward summation (alpha-pass) to a given sequence of 
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
            # prob[i,current_state] == P(O_i,q_i | lambda)
            for (past_state in state_label) {
                # P(q_i | q_{i-1}, lambda) * P(O_{i-1}, q_{i-1} | lambda)
                increment <- trans_mat[past_state,current_state] * probs[i-1,past_state]
                probs[i,current_state] <- probs[i,current_state] + increment
            }
            # Multiply by P(o_i | q_i, lambda)
            probs[i,current_state] <- probs[i,current_state] * outcome_mat[current_state,obs[i]]
        }
    }
    return(probs)
}

#------------------------------------------------------------------------------
# Description
#   Perform backward algorithm (beta-pass) to a given sequence of 
#   observations to probabilities of each state on all the stages
# Arguments:
#   obs: Sequence of observations
#   trans_mat: The transition matrix P
#   outcome_mat: The observation matrix E
#------------------------------------------------------------------------------
backward_sum <- function(obs, trans_mat, outcome_mat) {
    # Store some quantities for convenience
    len <- length(obs)    # Length of observations
    state_label <- colnames(trans_mat)   # The label of states, e.g., {"H","T"}
    num_state <- length(state_label)    # The number of states
    # Set up objects to store the results
    probs <- matrix(0, nrow=len, ncol=num_state)
    probs[len,] <- 1
    colnames(probs) <- state_label
    # Start calculating
    for (i in (len - 1):1) {
        for (current_state in state_label) {
            # probs[i,current_state] == P(o_{i+1},...,o_t | q_i, lambda)
            for (future_state in state_label) {
                # P(q_{i+1} | q_i, lambda) * P(o_{i+1} | q_{i+1}, lambda) * P(o_{i+2},...,o_t | q_{i+1}, lambda)
                increment <- trans_mat[current_state,future_state] * outcome_mat[future_state,obs[i+1]] * probs[i+1,future_state]
                probs[i,current_state] <- probs[i,current_state] + increment
            }
        }
    }
    return(probs)
}

#------------------------------------------------------------------------------
# Description
#   Combine forward summation (alpha-pass) and backward summation (beta-pass) 
#   to find an optimal sequence for a given sequence of observations
# Arguments:
#   obs: Sequence of observations
#   init_prob: The initial probability
#   trans_mat: The transition matrix P
#   outcome_mat: The observation matrix E
#------------------------------------------------------------------------------
forward_backward <- function(obs, init_prob, trans_mat, outcome_mat) {
    len <- length(obs)
    prob_alpha <- forward_sum(obs, init_prob, trans_mat, outcome_mat)   # The probability calculated by alpha-pass
    prob_beta <- backward_sum(obs, trans_mat, outcome_mat)   # The probability calculated by alpha-pass
    prob_t <- sum(prob_alpha[len,])     # The cumulative probability
    likelihood <- prob_alpha * prob_beta / prob_t     # The likelihood for each stage
    return(colnames(likelihood)[apply(likelihood, 1, which.max)])
}

#------------------------------------------------------------------------------
# Description
#   Use Viterbi algorithm to find an optimal state sequence for a given 
#   sequence of observations
# Arguments:
#   obs: Sequence of observations
#   init_prob: The initial probability
#   trans_mat: The transition matrix P
#   outcome_mat: The observation matrix E
#------------------------------------------------------------------------------
viterbi <- function(obs, init_prob, trans_mat, outcome_mat) {
    # Store some quantities for convenience
    len <- length(obs)
    state_label <- colnames(trans_mat)   # The label of states, e.g., {"H","T"}
    num_state <- length(state_label)    # The number of states
    # Initialize the matrices v and s
    v <- matrix(0, nrow=len, ncol=num_state)
    s <- matrix(0, nrow=len, ncol=num_state)
    colnames(v) <- state_label
    colnames(s) <- state_label
    # Fill the first row
    for (k in 1:num_state) {
        v[1, k] <- log(outcome_mat[k, obs[1]] * init_prob[k])
    }
    # Continue filling
    for (i in 2:len) {
        for (current_state in state_label) {
            tmp_prob <- matrix(0, nrow=1, ncol=num_state)
            colnames(tmp_prob) <- state_label
            for (past_state in state_label) {
                tmp_prob[1, past_state] <- log(trans_mat[past_state, current_state]) + v[i-1, past_state]
            }
            v[i, current_state] <- log(outcome_mat[current_state, obs[i]]) + max(tmp_prob)
            s[i, current_state] <- state_label[which.max(tmp_prob)]
        }
    }
    # Back tracing the states
    Q <- c(state_label[which.max(v[len,])])
    for (i in len:2) {
        Q <- c(as.character(s[len-1, Q[1]]), Q)
    }
    return(Q)
}
