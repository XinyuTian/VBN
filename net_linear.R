VBML_net <- function(x, y, L, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  if(is.null(parameters$c0)) parameters$c0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  
  N = nrow(x)
  P = ncol(x)
  
  alpha = list(a = parameters$a0 + P / 2, b = parameters$b0)
  lambda = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  
  # <beta * beta>
  E_betaSq = sum(beta$mu ^ 2) + sum(diag(beta$var))
  # <lambda>
  E_lambda = lambda$c / lambda$d
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <beta>
  E_beta = beta$mu
  # X* X
  xSq = crossprod(x, x)
  
  if(rec) {
    alpha.rec <- matrix(unlist(alpha), nrow = 1)
    beta.rec <- matrix(beta$mu, nrow = 1)
    lambda.rec <- matrix(unlist(lambda), nrow = 1)
  }
  
  for (i in 1 : (niter-1)) {
    # update alpha
    alpha$b = parameters$b0 + E_betaSq / 2
    E_alpha = drop(alpha$a / alpha$b)
    
    # update beta
    beta$var = diag(E_alpha, P, P) + E_lambda * xSq
    beta$mu = E_lambda * solve(beta$var) %*% t(x) %*% y
    E_betaSq = sum(beta$mu ^ 2) + sum(diag(beta$var))
    E_beta = beta$mu
    
    # update lambda
    lambda$d = parameters$d0 + crossprod(y, y) / 2 - t(E_beta) %*% t(x) %*% y + 0.5 * t(E_beta) %*% t(x) %*% x %*% E_beta
    E_lambda = drop(lambda$c / lambda$d)
    
    if(rec) {
      alpha.rec <- rbind(alpha.rec, unlist(alpha))
      beta.rec <- rbind(beta.rec, t(beta$mu)) 
      lambda.rec <- rbind(lambda.rec, unlist(lambda))
    }
  }
  
  if(rec){
    return(list(alpha = alpha, beta = beta, lambda = lambda, alpha.rec = alpha.rec, beta.rec = beta.rec, lambda.rec = lambda.rec))
  }
  return(list(alpha = alpha, beta = beta, lambda = lambda))
}