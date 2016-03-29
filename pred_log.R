pred_log = function(x, y, beta, npara = NULL) {
  # input y = {1, -1}
  if(!is.null(npara)) {
    beta0 = beta[1]
    beta = beta[-1]
    beta[(order(abs(beta), decreasing = T))[(npara+1):length(beta)]] = 0
    beta = matrix(c(beta0, beta), ncol = 1)
  }
  
  y = y / 2 + 0.5
  prob = sigmoid(cbind(1,x) %*% matrix(beta, ncol = 1))
  ypred = round(prob)
  bscore <- mean((y - prob)^2)
  acc <- sum(y == ypred) / length(y)
  return(list(brier = bscore, accuracy = acc))
}

comp_log <- function(x, y, x.test, y.test, Lmatrix=NULL) {
  ## ---------------------- MODEL ---------------- ##
  prior = list(a0 = runif(1, 0, 10), b0 = runif(1, 0, 10))
  res = glm((y / 2 + 0.5) ~ x, family = binomial)
  pred_log(x, y, res$coef)
  glm.fit = cv.glmnet(x, (y / 2 + 0.5), family = "binomial")
  beta.1se <- beta_glm(glm.fit, type = "1se") 
  pred_log(x, y, beta.1se)
  beta.min <- beta_glm(glm.fit, type = "min") 
  pred_log(x, y, beta.min)
  
  param = list(a0 = runif(1,0,10), b0 = runif(1,0,10))
  res = log_ridge(x, y, param)
  pred_log(x, y, res$beta$mu, npara = 9)
  
  res = VBML_ridge(x, y, prior, rec=T)
  
  glm.fit = cv.glmnet(x, y)
  
  res_net <- VBML_net(x, y, Lmatrix)
  res_lasso <- VBML_lasso(x, y, rec=T)
  res_lasso2 <- VBML_lasso2(x, y, rec=T)
  ## ---------------------- PREDICTION ---------------- ##
  beta1 <- beta_glm(glm.fit, type = "1se") 
  ssr1 <- ssr.fn(beta1, x.test, y.test)
  k1 = sum(beta1!=0)
  ssr0 <- ssr.fn(res$beta$mu, x.test, y.test, npara = k1)
  ssrN <- ssr.fn(res_net$beta$mu, x.test, y.test, npara = k1)
  ssrL <- ssr.fn(res_lasso$beta$mu, x.test, y.test, npara = k1)
  ssrL2 <- ssr.fn(res_lasso2$beta$mu, x.test, y.test, npara = NULL)
  w <- sqrt(2*pi*(res_lasso2$alpha$a / res_lasso2$alpha$b)/(res_lasso2$delta$c / res_lasso2$delta$d)) / 10
  comp_1se <- c(ridge=ssr0, Lasso=ssrL, network=ssrN, Lasso2=ssrL2, glmnet=ssr1, k=k1, kl=sum(res_lasso2$beta$mu!=0)-1, w=w)
  
  beta2 <- beta_glm(glm.fit, type = "min")
  ssr2 <- ssr.fn(beta2, x.test, y.test)
  k2 = sum(beta2!=0)
  ssr0 <- ssr.fn(res$beta$mu, x.test, y.test, npara = k2)
  ssrN <- ssr.fn(res_net$beta$mu, x.test, y.test, npara = k2)
  ssrL <- ssr.fn(res_lasso$beta$mu, x.test, y.test, npara = k2)
  ssrL2 <- ssr.fn(res_lasso2$beta$mu, x.test, y.test, npara = NULL)
  comp_min <- c(ridge=ssr0, Lasso=ssrL, network=ssrN, Lasso2=ssrL2, glmnet=ssr2, k=k2, kl=sum(res_lasso2$beta$mu!=0)-1, w=w)
  
  return(list('1se'=comp_1se, 'min'=comp_min))
}
