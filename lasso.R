# based on The Bayesian Lasso by T. Park and G. Casella
# gibbs -> variational Bayesian

VBML_lasso <- function(x, y, parameters = NULL, niter = 100, seed = NULL) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  if(is.null(parameters$c0)) parameters$c0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  if(is.null(parameters$e0)) parameters$e0 = 1
  if(is.null(parameters$f0)) parameters$f0 = 1
  
  if(! is.null(parameters$lambda0)) {
    if (is.numeric(parameters$lambda0)) {
      parameters$lambda.method = 'NU'
      warning("lambda would not be updated")
    } else {
      parameters$lambda0 = NULL
      warning("non-numeric lambda0 is dismissed")
    }
  }
  
  if(is.null(parameters$lambda.method)) {
    parameters$lambda.method = 'NU'
  } else if(! parameters$lambda.method %in% c("hyperprior", "empirical", "NU")) {
    parameters$lambda.method = "hyperprior"
    warning("lambda will be updated using a 'hyperprior'")
  }
  
  if(is.null(parameters$lambda0)) {
    parameters$lambda0 = switch(parameters$lambda.method, 
                                "NU" = 1, 
                                'empirical' = init.lambda(x, y),
                                'hyperprior' = parameters$e0 / parameters$f0
                                )
  }
  N = nrow(x)
  P = ncol(x)
  
  alpha = list(a = parameters$a0 + P / 2, b = parameters$b0)
  delta = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  gamma = list(lambda = parameters$lambda0)
  if (parameters$lambda.method == "hyperprior") {
    lambda = list(e = parameters$e0 + P, f = parameters$f0)
  }
  # <beta_j ^ 2>
  E_betajSq = beta$mu ^ 2 + diag(beta$var)
  # <beta * beta>
  E_betaSq = sum(E_betajSq)
  # <delta>
  E_delta = delta$c / delta$d
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <beta>
  E_beta = beta$mu
  # <gamma_j>
  E_gammaj = sqrt(gamma$lambda / E_alpha / E_betajSq)
  # <lambda>
  if (parameters$lambda.method == "hyperprior") {
    E_lambda = gamma$lambda = lambda$e / lambda$f
  }
  
  # X* X
  xSq = crossprod(x, x)
  
  for (i in 1 : (niter-1)) {
    # update alpha
    alpha$b = parameters$b0 + 0.5 * sum(E_gammaj * E_betajSq)
    E_alpha = drop(alpha$a / alpha$b)
    
    # update lambda
    if (parameters$lambda.method == "hyperprior") {
      lambda$f = parameters$f0 + sum(1 / E_gammaj) / 2
      E_lambda = gamma$lambda = lambda$e / lambda$f
    }
    # update gamma
    E_gammaj = sqrt(gamma$lambda / E_alpha / E_betajSq)
    
    # update beta
    beta$var = diag(drop(E_alpha * E_gammaj), P, P) + E_delta * xSq
    beta$mu = E_delta * solve(beta$var) %*% t(x) %*% y
    E_betajSq = beta$mu ^ 2 + diag(beta$var)
    E_betaSq = sum(E_betajSq)
    E_beta = beta$mu
    
    # update delta
    delta$d = parameters$d0 + crossprod(y, y) / 2 - t(E_beta) %*% t(x) %*% y + 0.5 * t(E_beta) %*% t(x) %*% x %*% E_beta
    E_delta = drop(delta$c / delta$d)
    
    # if(i < 5) print(beta$mu)
  }
  if (parameters$lambda.method == "hyperprior") {
    out = list(alpha = alpha, gamma = gamma, beta = beta, delta = delta, lambda = lambda)
  } else out = list(alpha = alpha, gamma = gamma, beta = beta, delta = delta)
  return(out)
}


init.lambda <- function(x, y) {
  lr <- lm(y ~ -1 + x)
  return( P * summary(lr)$sigma / sum(abs(lr$coefficients)))
}