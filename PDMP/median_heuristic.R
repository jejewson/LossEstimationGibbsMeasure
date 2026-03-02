median_heuristic <- function(Z) {
  # Compute the Gaussian kernel hyperparameter by median heuristic
  
  size1 <- length(Z)
  
  if (size1 > 100) {
    Zmed <- Z[1:100, , drop = FALSE]
    size1 <- 100
  } else {
    Zmed <- Z
  }
  
  # Compute squared norms of each row
  G <- rowSums(Zmed * Zmed)
  
  # Create matrices Q and R by replicating G
  Q <- matrix(rep(G, size1), nrow = size1, byrow = FALSE)
  R <- matrix(rep(G, each = size1), nrow = size1, byrow = TRUE)
  
  # Compute pairwise squared distances
  dists <- Q + R - 2 * (Zmed %*% t(Zmed))
  
  # Zero out lower triangle (including diagonal)
  dists <- dists - lower.tri(dists, diag = TRUE) * dists
  
  # Reshape to vector
  dists <- as.vector(dists)
  
  # Compute sigma as sqrt of half the median of positive distances
  sigma <- sqrt(0.5 * median(dists[dists > 0]))
  
  return(sigma)
}
