# Calculating Hotelling T-Square Statistic for a Test of Mean Vector
# arguments: data -- data set
#            mu0 -- mean vector in H0
#            alpha -- significant level
HotellingTSq <- function(data, mu0, alpha) {
  n <- nrow(data)
  p <- ncol(data)
  xbar <- apply(data, 2, mean)
  S <- var(data)
  Tsq <- n * t(xbar - mu0) %*% solve(S) %*% (xbar - mu0)
  T.critical <- ((n - 1) * p / (n - p)) * qf(1 - alpha, p, n - p)
  # outputs
  cat("n =", n, ", p =", p, "\n")
  cat("The T-square statistic is", round(Tsq, 3), "\n")
  cat("The critical value is", round(T.critical, 3), "\n")
  if(Tsq > T.critical) {
    cat("We reject H0 at significant level", 100 * (1 - alpha), "%")
  } else {
    cat("We fail to reject H0 at significant level", 100 * (1 - alpha), "%")
  }
}

# Obtaining Simultaneous confidence intervals for a mean vector
# arguments: data -- data set
#            alpha -- significant level
#            options -- types of confidence intervals
MultiIntervals <- function(data, alpha, options = "tsquare") {
  if (length(alpha) > 1) {
    stop("The significant level must be a number between 0 and 1!")
  }
  n <- nrow(data)
  p <- ncol(data)
  xbar <- apply(data, 2, mean)
  S <- var(data)
  if (options == "tsquare") {
    critical <- ((n - 1) * p / (n - p)) * qf(1 - alpha, p, n - p)
    cat("The simultaneous T-square", 100 * (1 - alpha), "% C.I.'s are:\n")
  } else if (options == "largesample") {
    critical <- qchisq(1 - alpha, p)
    cat("The large sample", 100 * (1 - alpha), "% C.I.'s are:\n")
  } else if (options == "bonferroni") {
    critical <- qt(1 - alpha / (2 * p), n - 1)^2
    cat("The Bonferroni", 100 * (1 - alpha), "% C.I.'s are:\n")
  } else {
    stop("The options must be one of the followings: tsquare, largesample, bonferroni")
  }
  for (i in 1:p) {
    moe <- sqrt(critical * S[i, i] / n)
    cat(i, ": (", round(xbar[i] - moe, 3), ",", round(xbar[i] + moe, 3), ")\n")
  }
}

chisqplot <- function(data, title = "Chi-Square Plot") {
  xbar <- apply(data, 2, mean)
  S <- cov(data)
  n <- nrow(data)
  p <- ncol(data)
  diff <- data.matrix(data - rep(1, n) %*% t(xbar))
  d <- diag(diff %*% solve(S) %*% t(diff))
  q <- qchisq((1:n - 0.5) / n, p)
  plot(sort(d)~q, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = title)
  abline(a = 0, b = 1)
}

compareMean <- function(data1, data2, diffx, alpha, conf.int = FALSE, conf.option = "tsquare") {
  if (ncol(data1) != ncol(data2)) {
    stop("The numbers of variables in these two data sets do not match!")
  }
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)
  xbar1 <- apply(data1, 2, mean)
  xbar2 <- apply(data2, 2, mean)
  diff <- xbar1 - xbar2
  S1 <- var(data1)
  S2 <- var(data2)
  S.pooled <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  
  S <- S.pooled * (1 / n1 + 1 / n2)
  Tsq <- t(diff - diffx) %*% solve(S) %*% (diff - diffx)
  c <- (n1 + n2 - 2) * p * qf(1 - alpha, p, n1 + n2 - p - 1) / (n1 + n2 - p - 1)
  
  cat("n1 =", n1, ", "n2 =", n2, ", p =", p, "\n")
  cat("The T-square statistic is", round(Tsq, 3), "\n")
  cat("The critical value is", round(c, 3), "\n")
  if(Tsq > c) {
    cat("We reject H0 at significant level", 100 * (1 - alpha), "%")
  } else {
    cat("We fail to reject H0 at significant level", 100 * (1 - alpha), "%")
  }
  
  if (!conf.int) {
    if (conf.option == "tsquare") {
      cat("The simultaneous T-square", 100 * (1 - alpha), "% C.I.'s are:\n")
      for (i in 1:p) {
        moe <- sqrt(c * S[i, i])
        cat(i, ": (", round(diff[i] - moe, 3), ",", round(diff[i] + moe, 3), ")\n")
      }
    } else if (conf.option == "bonferroni") {
      cat("The Bonferroni", 100 * (1 - alpha), "% C.I.'s are:\n")
      t <- qt(1 - alpha / (2 * p), n1 + n2 - 2)
      for (i in 1:p) {
        moe <- t * sqrt(S[i, i])
        cat(i, ": (", round(diff[i] - moe, 3), ",", round(diff[i] + moe, 3), ")\n")
      }
    } else {
      stop("The options must be one of the followings: tsquare, bonferroni")
    }
  }
}

normplot <- function(v, title) {
  n <- length(v)
  q <- qnorm((1:n - 0.5) / n)
  r_Q <- cor(q, v)
  
  qqnorm(v, main = title)
  qqline(v)
  
  r <- format(r_Q, digits = 3)
  eq <- bquote(r[Q] == .(r))
  mtext(eq)
}

ldiscriminant <- function(data1, data2) {
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  xbar1 <- apply(data1, 2, mean)
  xbar2 <- apply(data2, 2, mean)
  S1 <- var(data1)
  S2 <- var(data2)
  Spooled <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  a <- t(xbar1 - xbar2) %*% solve(Spooled)
  m <- as.numeric(a %*% (xbar1 + xbar2) / 2)
  
  confusion <- matrix(rep(0, 4), nrow = 2, ncol = 2)
  rownames(confusion) <- c("population 1", "population 2")
  colnames(confusion) <- c("population 1", "population 2")
  confusion[1, 1] <- sum(a %*% t(data.matrix(data1)) >= m)
  confusion[2, 2] <- sum(a %*% t(data.matrix(data2)) < m)
  confusion[1, 2] <- n1 - confusion[1, 1]
  confusion[2, 1] <- n2 - confusion[2, 2]
  return(confusion)
}

classification <- function(data1, data2, obs) {
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  xbar1 <- apply(data1, 2, mean)
  xbar2 <- apply(data2, 2, mean)
  Spooled <- ((n1 - 1) * var(data1) + (n2 - 1) * var(data2)) / (n1 + n2 - 2)
  a <- t(xbar1 - xbar2) %*% solve(Spooled)
  m <- as.numeric(a %*% (xbar1 + xbar2) / 2)
  if (as.numeric(a %*% t(obs)) >= m) {
    return(1)
  } else {
    return(2)
  }
}


holdout <- function(data1, data2) {
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  pop1 <- rep(0, n1)
  for (i in 1:n1) {
    pop1[i] <- classification(data1[-i,], data2, data1[i,])
  }
  pop2 <- rep(0, n2)
  for (j in 1:n2) {
    pop2[j] <- classification(data1, data2[-j,], data2[j,])
  }
  
  confusion <- matrix(rep(0, 4), nrow = 2, ncol = 2)
  rownames(confusion) <- c("population 1", "population 2")
  colnames(confusion) <- c("population 1", "population 2")
  confusion[1, 1] <- sum(pop1 == 1)
  confusion[1, 2] <- sum(pop1 == 2)
  confusion[2, 1] <- sum(pop2 == 1)
  confusion[2, 2] <- sum(pop2 == 2)

  return(confusion)
}


# Output the summary of clusters
# arguments: clusters -- output from cutree
clustersummary <- function(clusters) {
  M <- matrix(rep(0, 4*4), nrow = 4, ncol = 4)
  rownames(M) <- c("Cluster 1:", "Cluster 2:", "Cluster 3:", "Cluster 4:")
  colnames(M) <- c("1", "3", "6", "10")
  for (i in 1:4) {
    distr <- summary(factor(as.numeric(names(which(clusters == i))) %% 11))
    for (j in c("1", "3", "6", "10")) {
      M[i, j] <- ifelse(is.na(distr[j]), 0, distr[j])
    }
  }
  print(M)
}

clustersummary <- function(clusters) {
    labels <- c("1", "3", "6", "10")
    M <- matrix(rep(0, 4*4), nrow = 4, ncol = 4)
    rownames(M) <- c("Cluster 1:", "Cluster 2:", "Cluster 3:", "Cluster 4:")
    colnames(M) <- labels
    for (i in 1:4) {
        # We use the row number to determine the actual class of each observation in the cluster
        distr <- summary(factor(as.numeric(names(which(clusters == i))) %% 11))
        for (j in labels) 
            M[i, j] <- ifelse(is.na(distr[j]), 0, distr[j])
    }
    # Permute the labels and find the accuracy
    labels_perm <- permn(labels)
    max_corr <- 0
    max_index <- 0
    for (n in 1:factorial(4)) {
        rownames(M) <- labels_perm[[n]]
        correction <- sum(M[matrix(rep(labels, 2), nrow = 4, ncol = 2)])
        max_index <- ifelse(max_corr >= correction, max_index, n)
        max_corr <- ifelse(max_corr >= correction, max_corr, correction)
    }
    cat("Accuracy =", max_corr / 192)
    cat(", Classification:\n")
    rownames(M) <- labels_perm[[max_index]]
    M[1:4,] <- M[labels,]
    print(M)
}
