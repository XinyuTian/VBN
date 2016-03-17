# based on The Bayesian Lasso by T. Park and G. Casella
# gibbs -> variational Bayesian

VBML_lasso <- function(x, y, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  if(is.null(parameters$c0)) parameters$c0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  x = cbind(1, x)
  
#   if(! is.null(parameters$lambda0)) {
#     if (is.numeric(parameters$lambda0)) {
#       parameters$lambda.method = 'NU'
#       warning("lambda would not be updated")
#     } else {
#       parameters$lambda0 = NULL
#       warning("non-numeric lambda0 is dismissed")
#     }
#   }
#   
#   if(is.null(parameters$lambda.method)) {
#     parameters$lambda.method = 'NU'
#   } else if(! parameters$lambda.method %in% c("hyperprior", "empirical", "NU")) {
#     parameters$lambda.method = "hyperprior"
#     warning("lambda will be updated using a 'hyperprior'")
#   }
#   
#   if(is.null(parameters$lambda0)) {
#     parameters$lambda0 = switch(parameters$lambda.method, 
#                                 "NU" = 1, 
#                                 'empirical' = init.lambda(x, y),
#                                 'hyperprior' = parameters$e0 / parameters$f0
#                                 )
#   }
  N = nrow(x)
  P = ncol(x)
  
  alpha = list(a = parameters$a0 + P, b = parameters$b0)
  delta = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  gamma = list(mu = matrix(sqrt(alpha$a / alpha$b), nrow = P, ncol = 1), lambda = alpha$a / alpha$b)
#   if (parameters$lambda.method == "hyperprior") {
#     lambda = list(e = parameters$e0 + P, f = parameters$f0)
#   }
  # <beta_j ^ 2>
  E_betajSq = beta$mu ^ 2 + diag(beta$var)
  # <beta>
  E_beta = beta$mu
  # <delta>
  E_delta = delta$c / delta$d
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <gamma_j>
  E_gammaj = gamma$mu
  # <tau_j ^ 2>
  E_tauj = 1 / gamma$mu + 1/ gamma$lambda
  
  # X* X
  xSq = crossprod(x, x)
  
  if(rec) {
    alpha.rec <- matrix(unlist(alpha), nrow = 1)
    beta.rec <- matrix(beta$mu, nrow = 1)
    delta.rec <- matrix(unlist(delta), nrow = 1)
    gamma.rec <- matrix(unlist(gamma), nrow = 1)
  }

  for (i in 1 : (niter-1)) {
    # update alpha
    alpha$b = parameters$b0 + 0.5 * sum(E_tauj)
    E_alpha = drop(alpha$a / alpha$b)
    
    # update delta
    delta$d = parameters$d0 + crossprod(y, y) / 2 - t(E_beta) %*% t(x) %*% y + 0.5 * t(E_beta) %*% xSq %*% E_beta
    E_delta = drop(delta$c / delta$d)
    
    # update gamma
    gamma$mu = sqrt(E_alpha / E_betajSq)
    gamma$lambda = E_alpha
    E_gammaj = gamma$mu
    E_tauj = 1 / gamma$mu + 1/ gamma$lambda
    
    # update beta
    beta$delta = diag(E_gammaj, P, P) + E_delta * xSq
    beta$var = solve(beta$delta)
    beta$mu = E_delta * beta$var %*% t(x) %*% y
    E_betajSq = beta$mu ^ 2 + diag(beta$var)
    E_beta = beta$mu
    
    if(rec) {
      alpha.rec <- rbind(alpha.rec, unlist(alpha))
      beta.rec <- rbind(beta.rec, t(beta$mu)) 
      delta.rec <- rbind(delta.rec, unlist(delta))
      gamma.rec <- rbind(gamma.rec, unlist(gamma))
      
    }
  }
  if(rec){
    return(list(alpha = alpha, beta = beta, gamma = gamma, delta = delta, alpha.rec = alpha.rec, beta.rec = beta.rec, gamma.rec = gamma.rec, delta.rec = delta.rec))
  }
  return(list(alpha = alpha, beta = beta, gamma = gamma, delta = delta))
}


# init.lambda <- function(x, y) {
#   lr <- lm(y ~ -1 + x)
#   return( P * summary(lr)$sigma / sum(abs(lr$coefficients)))
# }