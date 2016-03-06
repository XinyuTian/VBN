## variational Bayesian multiple linear regression
rm(list=ls())
library(MASS)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_linear.R')
source('ridge.R')

## --------------------- SIMULATION -------------- ##
N = 100 # number of samples
P = 30 # number of features
aP = 5 # number of active features

a = 2; b = 1
c = 1; d = 4

data = simulation(N, P, aP, a, b, c, d, eff.size = 2)
x = data[[1]]
y = data[[2]]

# test with glmnet
# res = cv.glmnet(x, y)
## -------------------------  END  ------------- ##


## ---------------------- MODEL ---------------- ##
parameters = list(a0 = runif(1, 0, 100), b0 = runif(1, 0, 100), c0 = runif(1, 0, 100), d0 = runif(1, 0, 100))
res = VBML_ridge(x, y, parameters, rec=T)
#res$alpha
#res$lambda
#res$beta$mu

#res$alpha.rec

glm.fit = cv.glmnet(x, y)


## ---------------------- PREDICTION ---------------- ##
## prediction accuracy
acc.fn <- function(beta, x, y, npara = NULL) {
  if(!is.null(npara)) beta[order(abs(beta), decreasing = T)][(npara+1):length(beta)] = 0
    
  pred <- x %*%  beta
  acc <- sqrt(sum((pred - y)^2))
  return(acc)
}

beta_glm <- function(cv_fit, type = "1se") {
  lambda <- cv_fit[[paste0("lambda.", type)]]
  ind <- which(cv_fit$lambda == lambda)
  beta <- cv_fit$glmnet.fit$beta[, ind]
  return(beta)
}

acc0 <- acc.fn(res$beta$mu, x, y, npara = 7)

beta1 <- beta_glm(glm.fit) 
acc1 <- acc.fn(beta1, x, y)

beta2 <- beta_glm(glm.fit, type = "min")
acc2 <- acc.fn(beta2, x, y)

acc0; acc1; acc2
