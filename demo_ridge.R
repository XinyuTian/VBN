## variational Bayesian multiple linear regression
rm(list=ls())
library(MASS)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_linear.R')
source('pred_linear.R')
source('ridge.R')
source('lasso.R')
source('net_linear.R')
source('/media/slowsmile/lizhen_WD/Xinyu/NGLasso/simulation.R')

## --------------------- SIMULATION -------------- ##
N = 100 # number of samples
P = 30 # number of features
aP = 5 # number of active features

a = 2; b = 1
c = 1; d = 4

Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=5,ncol=5), rep=6)
Dmatrix <- diag(rowSums(Amatrix))
Lmatrix <- Dmatrix -Amatrix
di <- 1/sqrt(diag(Lmatrix))
Lmatrix <- t(t(Lmatrix*di)*di)

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

res_net <- VBML_net(x, y, Lmatrix)
res_lasso <- VBML_lasso(x, y, rec=T)
## ---------------------- PREDICTION ---------------- ##
## prediction accuracy

beta1 <- beta_glm(glm.fit) 
ssr1 <- ssr.fn(beta1, x, y)

beta2 <- beta_glm(glm.fit, type = "min")
ssr2 <- ssr.fn(beta2, x, y)

k = sum(beta2!=0)
ssr0 <- ssr.fn(res$beta$mu, x, y, npara = k)
ssrN <- ssr.fn(res_net$beta$mu, x, y, npara = k)
ssrL <- ssr.fn(res_lasso$beta$mu, x, y, npara = k)

ssr0; ssrN; ssrL; ssr1; ssr2; k
