simulation <- function(N, P, aP, a, b, c, d, eff.size = 1, seed = NULL) {
  if(! is.null(seed)) set.seed(seed)
  X <- matrix(nrow = N, ncol = P)
  y <- matrix(nrow = N, ncol = 1)
  for (i in 1:N) {
    alpha = rgamma(1, a, b)
    lambda = rgamma(1, c, d)
    beta = rnorm(P, 0, sd = sqrt(1 / alpha))
    beta[1:aP] = beta[1:aP] + eff.size
    
    X[i, ] = rnorm(P)
    
    y[i, ] = rnorm(1, t(beta[1:aP]) %*% X[i, 1:aP], sqrt(1 / lambda))
  }
  return (list("x" = X, "y" = y))
}