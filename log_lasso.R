log_lasso = function(x, y, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  
}

log_ridge = function(x, y, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  x = cbind(1, x)
  N = nrow(x)
  P = ncol(x)
  alpha = list(a = parameters$a0 + (P-1) / 2, b = parameters$b0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  xi = rnorm(N)
  
  # <beta * beta>
  E_betaSq = crossprod(beta$mu) + sum(diag(beta$var))
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <beta>
  E_beta = beta$mu
  # X* X
  xSq = crossprod(x, x)
  
  if(rec) {
    alpha.rec <- matrix(unlist(alpha), nrow = 1)
    beta.rec <- matrix(beta$mu, nrow = 1)
    xi.rec <- matrix(xi, nrow = 1)
  }
  
  for (i in 1 : (niter-1)) {
    # update alpha
    alpha$b = parameters$b0 + E_betaSq / 2
    E_alpha = drop(alpha$a / alpha$b)
    
    # update beta
    beta$delta = diag(E_alpha, P, P) - 2 * t(x) %*% diag(lam(xi)) %*% x 
    beta$var = solve(beta$delta)
    beta$mu = beta$var %*% t(x) %*% y / 2
    E_betaSq = sum(beta$mu ^ 2) + sum(diag(beta$var))
    E_beta = beta$mu
    
    # update xi
    xi = sqrt(diag(x %*% (beta$var + beta$mu %*% t(beta$mu)) %*% t(x)))
    
    if(rec) {
      alpha.rec <- rbind(alpha.rec, unlist(alpha))
      beta.rec <- rbind(beta.rec, t(beta$mu)) 
      xi.rec <- rbind(xi.rec, xi)
    }
  }
  
  if(rec){
    return(list(alpha = alpha, beta = beta, xi = xi, alpha.rec = alpha.rec, beta.rec = beta.rec, xi.rec = xi.rec))
  }
  return(list(alpha = alpha, beta = beta, xi = xi))
}


log_net = function(x, y, L, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  
}

