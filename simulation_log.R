simulation_log <- function(N, P, aP = P, a, b, eff.size = 1, eff = "random", seed = NULL) {
  if(! is.null(seed)) set.seed(seed)
  X <- matrix(nrow = N, ncol = P)
  y <- matrix(nrow = N, ncol = 1)
  if (eff == "random") {
    for (i in 1:N) {
      alpha = rgamma(1, a, b)
      beta = rnorm(P, 0, sd = sqrt(1 / alpha))
      beta[1:aP] = beta[1:aP] + eff.size
      
      X[i, ] = rnorm(P)
      
      y[i, ] = round(sigmoid(t(beta) %*% X[i, ])) * 2 - 1
    }
  } else if (eff == "fixed") {
    alpha = rgamma(1, a, b)
    beta = rnorm(P, 0, sd = sqrt(1 / alpha))
    beta[1:aP] = beta[1:aP] + eff.size
    for (i in 1:N) {
      X[i, ] = rnorm(P)
    }
    y <- round(sigmoid(X %*% beta)) * 2 - 1
  }
  return (list("x" = X, "y" = y))
}

