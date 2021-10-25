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
