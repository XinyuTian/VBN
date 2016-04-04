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
  
  alpha = list(a = parameters$a0 + P-1, b = parameters$b0)
  delta = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  gamma = list(mu = matrix(sqrt(alpha$a / alpha$b), nrow = P-1, ncol = 1), lambda = alpha$a / alpha$b)
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
    gamma$mu = sqrt(E_alpha / E_betajSq[-1])
    gamma$lambda = E_alpha
    E_gammaj = gamma$mu
    E_tauj = 1 / gamma$mu + 1/ gamma$lambda
    
    # update beta
    beta$delta = diag(c(0,E_gammaj), P, P) + E_delta * xSq
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

VBML_lasso2 <- function(x, y, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  if(is.null(parameters$c0)) parameters$c0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  if(is.null(parameters$w)) parameters$w = 1
  N = nrow(x)
  x = cbind(1, x)
  P = ncol(x)
  gamma <- rep(TRUE, P-1);gammause <- c(TRUE,gamma)
  
  alpha = list(a = parameters$a0 + P-1, b = parameters$b0)
  delta = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  itau = list(mu = matrix(sqrt(alpha$a / alpha$b), nrow = P-1, ncol = 1), lambda = alpha$a / alpha$b)
  # <beta_j ^ 2>
  E_betajSq = beta$mu ^ 2 + diag(beta$var)
  # <beta>
  E_beta = beta$mu
  # <delta>
  E_delta = delta$c / delta$d
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <itau_j>
  E_itauj = itau$mu
  # <tau_j ^ 2>
  E_tauj = 1 / itau$mu + 1/ itau$lambda
  
  # X* X
  xSq = crossprod(x, x)
  
  if(rec) {
    alpha.rec <- matrix(unlist(alpha), nrow = 1)
    beta.rec <- matrix(beta$mu, nrow = 1)
    delta.rec <- matrix(unlist(delta), nrow = 1)
    itau.rec <- matrix(unlist(itau), nrow = 1)
    gamma.rec <- matrix(gamma, nrow = 1)
    w.rec <- matrix(parameters$w, nrow = 1)
  }
  
  for (i in 1 : (niter-1)) {
    # update beta
    if(sum(gammause) > 1) {
      beta$delta = diag(c(0,E_itauj[gamma])) + E_delta * xSq[gammause,gammause]
      beta$var[gammause,gammause] = solve(beta$delta)
    } else beta$var[gammause,gammause] <- matrix(1 / (E_delta * N))
    beta$mu[gammause] = E_delta * beta$var[gammause,gammause] %*% t(x[,gammause]) %*% y
    beta$mu[!gammause] = 0
    E_betajSq = beta$mu ^ 2 + diag(beta$var)
    E_beta = beta$mu
    
    # update alpha
    alpha$a = parameters$a0 + sum(gamma)
    alpha$b = parameters$b0 + 0.5 * sum(E_tauj[gamma])
    E_alpha = drop(alpha$a / alpha$b)
    
    # update delta
    delta$d = parameters$d0 + crossprod(y, y) / 2 - t(E_beta[gammause]) %*% t(x[,gammause]) %*% y + 0.5 * t(E_beta[gammause]) %*% xSq[gammause,gammause] %*% E_beta[gammause]
    E_delta = drop(delta$c / delta$d)
    
    # update itau
    itau$mu = sqrt(E_alpha / E_betajSq[-1])
    itau$lambda = E_alpha
    E_itauj = itau$mu
    E_tauj = 1 / itau$mu + 1/ itau$lambda
    
    # updata gamma
    lambda <- 2 * sqrt(E_alpha) / E_delta
    y_diff <-y - x[, gammause, drop=F] %*% beta$mu[gammause, , drop=F]
    #parameters$w <- sqrt(2*pi*E_alpha/E_delta) / 10
    for(j in seq(2,P)) {
      if(gamma[j-1]) {
        yj_fit <- x[,j, drop=F] %*% beta$mu[j]
        y_diff2 <- y_diff + yj_fit
        rr <- parameters$w * exp(- E_delta / 2 * (crossprod(y_diff) - crossprod(y_diff2) + lambda * abs(beta$mu[j])))
      } else{
        betaj <- crossprod(x[,j],y_diff) / crossprod(x[,j])
        yj_fit <- x[,j, drop=F] %*% betaj
        rr <- parameters$w * exp(- E_delta / 2 * (crossprod(y_diff-yj_fit) - crossprod(y_diff) + lambda * abs(betaj)))
      }
      gamma[j-1] <- (rr>1)
    }
    gammause <- c(TRUE,gamma)
    
    if(rec) {
      alpha.rec <- rbind(alpha.rec, unlist(alpha))
      beta.rec <- rbind(beta.rec, t(beta$mu)) 
      delta.rec <- rbind(delta.rec, unlist(delta))
      itau.rec <- rbind(itau.rec, unlist(itau))
      gamma.rec <- rbind(gamma.rec, gamma)
      w.rec <- rbind(w.rec, parameters$w)
    }
  }
  out <- list(alpha = alpha, beta = beta, itau = itau, delta = delta, gamma = gamma, w=parameters$w)
  if(rec){
    out <- c(out, list(alpha.rec = alpha.rec, beta.rec = beta.rec, delta.rec = delta.rec, itau.rec = itau.rec, gamma.rec=gamma.rec, w.rec=w.rec))
  }
  return(out)
}

VBML_lasso3 <- function(x, y, parameters = NULL, niter = 100, seed = NULL, rec = FALSE) {
  if(! is.null(seed)) set.seed(seed)
  if(is.null(parameters$a0)) parameters$a0 = 1
  if(is.null(parameters$b0)) parameters$b0 = 1
  if(is.null(parameters$c0)) parameters$c0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  if(is.null(parameters$d0)) parameters$d0 = 1
  if(is.null(parameters$w)) parameters$w = 1
  
  muy = mean(y)
  sdy = sd(y)
  y = (y - muy) / sdy
  N = nrow(x)
  x = cbind(1, x)
  P = ncol(x)
  gamma <- rep(TRUE, P-1);gammause <- c(TRUE,gamma)
  
  alpha = list(a = parameters$a0 + P-1, b = parameters$b0)
  delta = list(c = parameters$c0 + N / 2, d = parameters$d0)
  beta = list(mu = matrix(0, nrow = P, ncol = 1), var = diag(1, P, P))
  itau = list(mu = matrix(sqrt(alpha$a / alpha$b), nrow = P-1, ncol = 1), lambda = alpha$a / alpha$b)
  # <beta_j ^ 2>
  E_betajSq = beta$mu ^ 2 + diag(beta$var)
  # <beta>
  E_beta = beta$mu
  # <delta>
  E_delta = delta$c / delta$d
  # <alpha>
  E_alpha = alpha$a / alpha$b
  # <itau_j>
  E_itauj = itau$mu
  # <tau_j ^ 2>
  E_tauj = 1 / itau$mu + 1/ itau$lambda
  
  # X* X
  xSq = crossprod(x, x)
  
  if(rec) {
    alpha.rec <- matrix(unlist(alpha), nrow = 1)
    beta.rec <- matrix(beta$mu, nrow = 1)
    delta.rec <- matrix(unlist(delta), nrow = 1)
    itau.rec <- matrix(unlist(itau), nrow = 1)
    gamma.rec <- matrix(gamma, nrow = 1)
  }
  
  for (i in 1 : (niter-1)) {
    # update beta
    if(sum(gammause) > 1) {
      beta$delta = diag(c(0,E_itauj[gamma])) + E_delta * xSq[gammause,gammause]
      beta$var[gammause,gammause] = solve(beta$delta)
    } else beta$var[gammause,gammause] <- matrix(1 / (E_delta * N))
    beta$mu[gammause] = E_delta * beta$var[gammause,gammause] %*% t(x[,gammause]) %*% y
    beta$mu[!gammause] = 0
    E_betajSq = beta$mu ^ 2 + diag(beta$var)
    E_beta = beta$mu
    
    # update alpha
    alpha$a = parameters$a0 + sum(gamma)
    alpha$b = parameters$b0 + 0.5 * sum(E_tauj[gamma])
    E_alpha = drop(alpha$a / alpha$b)
    
    # update delta
    delta$d = parameters$d0 + crossprod(y, y) / 2 - t(E_beta[gammause]) %*% t(x[,gammause]) %*% y + 0.5 * t(E_beta[gammause]) %*% xSq[gammause,gammause] %*% E_beta[gammause]
    E_delta = drop(delta$c / delta$d)
    
    # update itau
    itau$mu = sqrt(E_alpha / E_betajSq[-1])
    itau$lambda = E_alpha
    E_itauj = itau$mu
    E_tauj = 1 / itau$mu + 1/ itau$lambda
    
    # updata gamma
    lambda <- 2 * sqrt(E_alpha) / E_delta * sdy
    #parameters$w <- sqrt(2*pi*E_alpha/E_delta) / 10
    logdev0 <- logdev(x[, gammause, drop=F], y, lambda)
    for(j in seq(2,P)) {
      gamma1 <- gammause; gamma1[j] <- !gammause[j]
      logdev1 <- logdev(x[, gamma1, drop=F], y, lambda)
      if(gamma[j-1]) {
        rr <- parameters$w * exp(- E_delta / 2 * (logdev0 - logdev1))
      } else{
        rr <- parameters$w * exp(- E_delta / 2 * (logdev1 - logdev0))
      }
#       if((rr>1) != gamma[j-1]) {
#         gamma[j-1] <- (rr>1)
#         gammause <- c(TRUE,gamma)
#         logdev0 <- logdev(x[, gammause, drop=F], y, lambda)
#       }
      gamma[j-1] <- (rr>1)
    }
    gammause <- c(TRUE,gamma)
    
    if(rec) {
      alpha.rec <- rbind(alpha.rec, unlist(alpha))
      beta.rec <- rbind(beta.rec, t(beta$mu)) 
      delta.rec <- rbind(delta.rec, unlist(delta))
      itau.rec <- rbind(itau.rec, unlist(itau))
      gamma.rec <- rbind(gamma.rec, gamma)
    }
  }
  beta$mu <- sdy * beta$mu; beta$mu[1] <- beta$mu[1] + muy
  beta$var <- (sdy ^ 2) * beta$var
  out <- list(alpha = alpha, beta = beta, itau = itau, delta = delta, gamma = gamma)
  if(rec){
    out <- c(out, list(alpha.rec = alpha.rec, beta.rec = beta.rec, delta.rec = delta.rec, itau.rec = itau.rec, gamma.rec=gamma.rec))
  }
  return(out)
}

## function return - 2*(loglike -loglike(Null)) + lambda * penalty
logdev <- function(x, y, lambda) {
  if(ncol(x) == 1) return(0)
  if(ncol(x) == 2) {
    # -2*loglike(Null)
    llN <- crossprod(y - mean(y))
    # -2*loglike
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    beta <- pmax(abs(b) - c(0, lambda / 2 / nrow(y)), 0) * sign(b)
    ll <- crossprod(y - x %*% matrix(beta))
    out = - ll + llN + lambda * abs(beta[2])
    return(out)
  }
  glmfit <- glmnet(x[,-1], y, family = 'gaussian', lambda = lambda / nrow(y) / 2)
  beta <- as.matrix(glmfit$beta)
  out = - glmfit$nulldev * glmfit$dev.ratio + lambda * sum(abs(beta))
  return(out)
}
