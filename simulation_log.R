simulation_log <- function(N, P, aP = P, a, b) {
  X <- matrix(nrow = N, ncol = P)
  y <- matrix(nrow = N, ncol = 1)
  for (i in 1:N) {
    alpha = rgamma(1, a, b)
    beta = rnorm(P, 0, sd = sqrt(1 / alpha))
    beta[1:aP] = beta[1:aP] + 2
    
    X[i, ] = rnorm(P)
    
    y[i, ] = round(sigmoid(t(beta[1:aP]) %*% X[i, 1:aP])) * 2 - 1
  }
  return (list("x" = X, "y" = y))
}

