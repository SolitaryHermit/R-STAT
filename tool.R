M <- matrix(c(0.6, 0.4, 0.02, 0.98), 2, 2, TRUE)
rownames(M) <- c("Present", "Absent")
colnames(M) <- c("P", "N")
prior <- matrix(c(0.001, 0.999))
rownames(prior) <- c("Present", "Absent")
ep <- function(s) {
  data_present <- prod(M["Present", s])
  data_absent <- prod(M["Absent", s])
  print(paste("P(data | Present) =", data_present))
  print(paste("P(data | Absent) =", data_absent))
  if (data_present > data_absent) {
    print("MLE is Present")
  } else {
    print("MLE is Absent")
  }
  p_present <- data_present * prior["Present",1]
  p_absent <- data_absent * prior["Absent",1]
  post_present <- p_present / (p_present + p_absent)
  post_absent <- p_absent / (p_present + p_absent)
  print(paste("P(Present | data) =", post_present))
  print(paste("P(Absent | data) =", post_absent))
  if (post_present > post_absent) {
    print("MAP is Present")
  } else {
    print("MAP is Absent")
  }
}
