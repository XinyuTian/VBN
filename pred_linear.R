## calculate prediction accuracy
ssr.fn <- function(beta, x, y, npara = NULL) {
  beta0 = beta[1]
  beta = beta[-1]
  if(!is.null(npara)) beta[(order(abs(beta), decreasing = T))[(npara+1):length(beta)]] = 0
  beta = matrix(c(beta0, beta), ncol = 1)
  
  pred <- cbind(1, x) %*%  beta
  ssr <- sqrt(sum((pred - y)^2))
  return(ssr)
}

beta_glm <- function(cv_fit, type = "1se") {
  lambda <- cv_fit[[paste0("lambda.", type)]]
  ind <- which(cv_fit$lambda == lambda)
  beta <- c(cv_fit$glmnet.fit$a0[ind], cv_fit$glmnet.fit$beta[, ind])
  return(beta)
}
