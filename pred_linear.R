## calculate prediction accuracy
acc.fn <- function(beta, x, y) {
  pred <- x %*%  beta
  acc <- sqrt(sum((pred - y)^2))
  return(acc)
}

## calculate beta from cv.glmnet
beta_glm <- function(cv_fit, type = "1se") {
  lambda <- cv_fit[[paste0("lambda.", type)]]
  ind <- which(cv_fit$lambda == lambda)
  beta <- cv_fit$glmnet.fit$beta[, ind]
  return(beta)
}

