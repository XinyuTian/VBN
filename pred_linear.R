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

comp <- function(x, y, x.test, y.test, Lmatrix) {
  ## ---------------------- MODEL ---------------- ##
  prior = list(a0 = runif(1, 0, 10), b0 = runif(1, 0, 10), c0 = runif(1, 0, 10), d0 = runif(1, 0, 10))
  res = VBML_ridge(x, y, prior, rec=T)
  
  glm.fit = cv.glmnet(x, y)
  
  res_net <- VBML_net(x, y, Lmatrix)
  res_lasso <- VBML_lasso(x, y, rec=T)
  res_lasso2 <- VBML_lasso2(x, y, rec=T)
  res_lasso3 <- VBML_lasso3(x, y, rec=T)
  ## ---------------------- PREDICTION ---------------- ##
  beta1 <- beta_glm(glm.fit, type = "1se") 
  ssr1 <- ssr.fn(beta1, x.test, y.test)
  k1 = sum(beta1!=0)
  ssr0 <- ssr.fn(res$beta$mu, x.test, y.test, npara = k1)
  ssrN <- ssr.fn(res_net$beta$mu, x.test, y.test, npara = k1)
  ssrL <- ssr.fn(res_lasso$beta$mu, x.test, y.test, npara = k1)
  ssrL2 <- ssr.fn(res_lasso2$beta$mu, x.test, y.test, npara = NULL)
  ssrL3 <- ssr.fn(res_lasso3$beta$mu, x.test, y.test, npara = NULL)
  #w <- sqrt(2*pi*(res_lasso2$alpha$a / res_lasso2$alpha$b)/(res_lasso2$delta$c / res_lasso2$delta$d)) / 2
  comp_1se <- c(ridge=ssr0, Lasso=ssrL, network=ssrN, Lasso2=ssrL2, Lasso3=ssrL3, glmnet=ssr1, k=k1, kl2=sum(res_lasso2$beta$mu!=0)-1, kl3=sum(res_lasso3$beta$mu!=0)-1)
  
  beta2 <- beta_glm(glm.fit, type = "min")
  ssr2 <- ssr.fn(beta2, x.test, y.test)
  k2 = sum(beta2!=0)
  ssr0 <- ssr.fn(res$beta$mu, x.test, y.test, npara = k2)
  ssrN <- ssr.fn(res_net$beta$mu, x.test, y.test, npara = k2)
  ssrL <- ssr.fn(res_lasso$beta$mu, x.test, y.test, npara = k2)
  ssrL2 <- ssr.fn(res_lasso2$beta$mu, x.test, y.test, npara = NULL)
  comp_min <- c(ridge=ssr0, Lasso=ssrL, network=ssrN, Lasso2=ssrL2, Lasso3=ssrL3, glmnet=ssr2, k=k2, kl=sum(res_lasso2$beta$mu!=0)-1, kl3=sum(res_lasso3$beta$mu!=0)-1)
  
  return(list('1se'=comp_1se, 'min'=comp_min))
}

## lower bound function
lbound <- function(x, y, beta){
  return(crossprod(y - x %*% beta))
}
