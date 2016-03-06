logistic <- function(x, y, parameters = NULL, niter = 1000, seed = NULL) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1

  N = nrow(x)
  P = ncol(x)
  
  alpha = list(a = parameters$a0 + P / 2, b = parameters$b0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))

  # <beta_j ^ 2>
  E_betajSq = beta$mu ^ 2 + diag(beta$var)
  # <beta * beta>
  E_betaSq = sum(E_betajSq)
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # xi and lam_xi
  lam_xi = rep( 1/8, N)
  # Sxy
  Sxy = t(x) %*% y / 2
  
  for (i in 1 : (niter-1)) {
    # update alpha
    alpha$b = parameters$b0 + 0.5 * E_betaSq
    E_alpha = drop(alpha$a / alpha$b)
    
    # update beta
    beta$var = E_alpha * diag(P) + 2 * t(x) %*% diag(lam_xi) %*% x
    V = solve(beta$var)
    beta$mu = V %*% Sxy
    E_betajSq = beta$mu ^ 2 + diag(beta$var)
    E_betaSq = sum(E_betajSq)
    
    # update xi
    V_xi = V + tcrossprod(beta$mu)
    xi = sqrt(apply(x, 1, function(xx) matrix(xx, nrow = 1) %*% V_xi %*% matrix(xx, ncol = 1)))
    lam_xi = lam(xi)
  }
  
  return(list(alpha = alpha, beta = beta))
}

lam <- function(xi) {
  return((sigmoid(xi) - 0.5) / 2 / xi)
}

sigmoid <- function(x) {
  return( 1 / (1 + exp(-x)))
}