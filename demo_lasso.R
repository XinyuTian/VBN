## variational Bayesian multiple linear regression
rm(list=ls())
library(MASS)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_linear.R')
source('lasso.R')

## --------------------- SIMULATION -------------- ##
N = 100 # number of samples
P = 30 # number of features
aP = 5 # number of active features

a = 100; b = 50
c = 20; d = 10

data = simulation(N, P, aP, a, b, c, d)
x = data[[1]]
y = data[[2]]

# test with glmnet
# res = cv.glmnet(x, y)
## -------------------------  END  ------------- ##


## ---------------------- MODEL ---------------- ##
parameters = list(lambda.method = 'hyperprior')
res = VBML_lasso(x, y, parameters)
res$beta$mu
res$gamma
