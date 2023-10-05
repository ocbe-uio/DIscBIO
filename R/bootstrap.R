bootstrap <- function(data, n_reps) {
  results <- numeric(n_reps)
  n <- nrow(data)
  i <- matrix(sample(n, n * n_reps, replace = TRUE), n_reps, n)
  freqs <- t(apply(i, 1, tabulate, ncol(i)))
  for (rep in seq_len(n_reps)) {
    samp <- freqs[rep, ]
    jac <- JS(data, samp)
    results[rep] <- mean(jac)
  }
  return(results)
}

JS <- function(data, indices) {
  d <- data[indices, ]
  jac <- suppressMessages(distance(t(d), method = "jaccard"))
  jac1 <- 1 - jac
  JSmean <- mean(jac1)
  return(JSmean)
}
