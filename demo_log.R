## variational Bayesian multiple logistic regression
rm(list=ls())
library(MASS)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_log.R')
source('logistic.R')
source('pred_log.R')
library(glmnet)
source('/media/slowsmile/lizhen_WD/Xinyu/NGLasso/simulation.R')

## --------------------- SIMULATION -------------- ##
N = 200 # number of samples
r = 1 # #{training samples} / #{test samples}
P = 30 # number of features
aP = 5 # number of active features

a = 100; b = 50

Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=aP,ncol=aP), rep=P/aP)
Dmatrix <- diag(rowSums(Amatrix))
Lmatrix <- Dmatrix -Amatrix
di <- 1/sqrt(diag(Lmatrix))
Lmatrix <- t(t(Lmatrix*di)*di)

data = simulation_log(N, P, aP, a, b)
x = data[[1]]
y = data[[2]]

# test with glmnet
res = glm((y / 2 + 0.5) ~ -1 + x, family = binomial())
pred_log(x, y, res$coef)
res = cv.glmnet(x, (y / 2 + 0.5), family = "binomial")
beta = glmnet(x, (y / 2 + 0.5), family = "binomial", lambda = res$lambda.1se)$beta
pred_log(x, y, beta)

## ---------------------------------------------- ##
n_sim = 100
res_1se <- matrix(nrow = 0, ncol = 7)
res_min <- matrix(nrow = 0, ncol = 7)
