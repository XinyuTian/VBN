## variational Bayesian multiple logistic regression
rm(list=ls())
library(MASS)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_log.R')
source('logistic.R')
source('pred_log.R')

## --------------------- SIMULATION -------------- ##
N = 100 # number of samples
P = 15 # number of features
aP = 5 # number of active features

a = 100; b = 50

data = simulation_log(N, P, aP, a, b)
x = data[[1]]
y = data[[2]]

# test with glmnet
res = glm((y / 2 + 0.5) ~ -1 + x, family = binomial())
pred_log(x, y, res$coef)
res = cv.glmnet(x, (y / 2 + 0.5), family = "binomial")
beta = glmnet(x, (y / 2 + 0.5), family = "binomial", lambda = res$lambda.1se)$beta
pred_log(x, y, beta)

## -------------------- FITTING ---------------- ##
res = logistic(x, y, niter= 100)
res$beta$mu
pred_log(x, y, res$beta$mu)
