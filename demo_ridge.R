## variational Bayesian multiple linear regression
rm(list=ls())
library(MASS)
library(glmnet)
setwd("/media/slowsmile/lizhen_WD/Xinyu/VBML/")
source('simulation_linear.R')
source('pred_linear.R')
source('ridge.R')
source('lasso.R')
source('net_linear.R')
source('/media/slowsmile/lizhen_WD/Xinyu/NGLasso/simulation.R')

## --------------------- SIMULATION -------------- ##
N = 200 # number of samples
r = 1 # #{training samples} / #{test samples}
P = 30 # number of features
aP = 5 # number of active features

a = 2; b = 1
c = 1; d = 4

Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=5,ncol=5), rep=6)
Dmatrix <- diag(rowSums(Amatrix))
Lmatrix <- Dmatrix -Amatrix
di <- 1/sqrt(diag(Lmatrix))
Lmatrix <- t(t(Lmatrix*di)*di)

sim_data = simulation(N, P, aP, a, b, c, d, eff.size = 2)
x = sim_data[[1]][1:round(N*r/(1+r)), , drop=F]
y = sim_data[[2]][1:round(N*r/(1+r)), , drop=F]
x.test = sim_data[[1]][(round(N*r/(1+r)) + 1):N, , drop=F]
y.test = sim_data[[2]][(round(N*r/(1+r)) + 1):N, , drop=F]

# test with glmnet
# res = cv.glmnet(x, y)
## -------------------------  END  ------------- ##
n_sim = 100
res_1se <- matrix(nrow = 0, ncol = 5)
res_min <- matrix(nrow = 0, ncol = 5)

for (i in seq(n_sim)){
  sim_data = simulation(N, P, aP, a, b, c, d, eff.size = 2)
  x = sim_data[[1]][1:round(N*r/(1+r)), , drop=F]
  y = sim_data[[2]][1:round(N*r/(1+r)), , drop=F]
  x.test = sim_data[[1]][(round(N*r/(1+r)) + 1):N, , drop=F]
  y.test = sim_data[[2]][(round(N*r/(1+r)) + 1):N, , drop=F]

  res <- comp(x, y, x.test, y.test, Lmatrix)
  res_1se <- rbind(res_1se, res$'1se')
  res_min <- rbind(res_min, res$min)
}

jpeg(filename = "~/Dropbox/My R Code/VBML/comp_min.jpeg")
boxplot(res_min[,1:4], main = "Sums of Squares of Residules (min)")
dev.off()
jpeg(filename = "~/Dropbox/My R Code/VBML/comp_1se.jpeg")
boxplot(res_1se[,1:4], main = "Sums of Squares of Residules (1se)")
dev.off()



apply(res_min[,1:4], 2, tapply, INDEX=res_min[,5], FUN=mean)
