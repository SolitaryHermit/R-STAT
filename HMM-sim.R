#------------------------------------------------------------------------------
# Description
#   Simulation of Hidden Markov Model
# Arguments:
#   n: Length of observations
#   init_prob: The initial probability
#   trans_mat: The transition matrix P
#   outcome_mat: The emission matrix E
#------------------------------------------------------------------------------
hidden_markov <- function(T, init_dist, trans_mat, emit_mat) {
    # Record the constants
    num_states <- nrow(emit_mat)
    num_outcomes <- ncol(emit_mat)
    # Save the labels for states
    if (is.null(rownames(emit_mat))) {
        state_labels <- 1:num_states
    } else {
        state_labels <- rownames(emit_mat)
    }
    # Save the labels for observations
    if (is.null(colnames(emit_mat))) {
        obs_labels <- 1:num_outcomes
    } else {
        obs_labels <- colnames(emit_mat)
    }
    # Set up the vectors recording the states and observations
    states <- rep(0, T)
    observations <- rep(0, T)
    # Initialize by initial distribution and transition matrix
    states[1] <- sample(1:num_states, 1, prob=init_dist)
    observations[1] <- sample(1:num_outcomes, 1, prob=emit_mat[states[1],])
    # Complete the following simulation
    for (t in 2:T) {
        states[t] <- sample(1:num_states, 1, prob=trans_mat[states[t-1],])
        observations[t] <- sample(1:num_outcomes, 1, prob=emit_mat[states[t],])
    }
    # Set up the labels
    states <- state_labels[states]
    observations <- obs_labels[observations]
    return(list(Q=states, O=observations))
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
forward_algorithm <- function(observations, init_dist, trans_mat, emit_mat, output=NA) {
    # Record the constants
    T <- length(observations)
    num_states <- nrow(emit_mat)
    # Set up forward table
    alpha <- matrix(0, nrow=num_states, ncol=T)
    # Initialization
    alpha[, 1] <- emit_mat[, observations[1]] * init_dist
    # Recursion
    for (t in 2:T) {
        alpha[, t] <- emit_mat[, observations[t]] * (t(trans_mat) %*% alpha[, t-1])
    }
    # termination
    prob <- sum(alpha[, T])
    # Set up return values
    if (output == "p") {
        return(prob)
    } else if (output == "t") {
        return(alpha)
    }
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
backward_algorithm <- function(observations, init_dist, trans_mat, emit_mat, output=NA) {
    # Record the constants
    T <- length(observations)
    num_states <- nrow(emit_mat)
    # Set up forward table
    beta <- matrix(0, nrow=num_states, ncol=T)
    # Initialization
    beta[, T] <- 1
    # Recursion
    for (t in (T-1):1) {
        beta[, t] <- trans_mat %*% (emit_mat[, observations[t+1]] * beta[, t+1])
    }
    # termination
    prob <- init_dist %*% (emit_mat[, observations[1]] * beta[, 1])
    # Set up return values
    if (output == "p") {
        return(prob)
    } else if (output == "t") {
        return(beta)
    }
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
fb_decoding <- function(observations, init_dist, trans_mat, emit_mat, output=NA) {
    # Record the constants
    T <- length(observations)
    # Save the labels for states
    if (is.null(rownames(emit_mat))) {
        state_labels <- 1:num_states
    } else {
        state_labels <- rownames(emit_mat)
    }
    # Obtain forward and backward tables and total likelihood
    alpha <- forward_algorithm(observations, init_dist, trans_mat, emit_mat, output="t")
    beta <- backward_algorithm(observations, init_dist, trans_mat, emit_mat, output="t")
    total_prob <- sum(alpha[, T])
    # total_prob <- init_dist %*% (emit_mat[, observations[1]] * beta[, 1])
    # Compute likelihood for each state at time t
    gamma <- alpha * beta / total_prob
    # Set up return values
    if (output == "q") {
        return(state_labels[apply(gamma, 2, which.max)])
    } else if (output == "t") {
        return(gamma)
    }
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
viterbi <- function(observations, init_dist, trans_mat, emit_mat, output=NA) {
    # Record the constants
    T <- length(observations)
    num_states <- nrow(emit_mat)
    # Save the labels for states
    if (is.null(rownames(emit_mat))) {
        state_labels <- 1:num_states
    } else {
        state_labels <- rownames(emit_mat)
    }
    # Set up matrix V and S
    V <- matrix(0, nrow=num_states, ncol=T)
    S <- matrix(0, nrow=num_states, ncol=T)
    # Initialization
    V[, 1] <- log(init_dist * emit_mat[, observations[1]])
    # Recursion
    for (t in 2:T) {
        for (s in 1:num_states) {
            values <- log(trans_mat[, s]) + V[, t-1]
            S[s, t] <- which.max(values)
            V[s, t] <- log(emit_mat[s, observations[t]]) + values[S[s, t]]
        }
    }
    # Termination
    states <- c(which.max[V[, T]])
    for (t in T:2) {
        states <- c(V[states[0], t], states)
    }
    states <- state_labels[states]
    # Set up return values
    if (output == "q") {
        return(states)
    } else if (output == "v") {
        return(V)
    } else if (output == "s") {
        return(S)
    }
}
