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

Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=aP,ncol=aP), rep=P/aP)
Dmatrix <- diag(rowSums(Amatrix))
Lmatrix <- Dmatrix -Amatrix
di <- 1/sqrt(diag(Lmatrix))
Lmatrix <- t(t(Lmatrix*di)*di)

sim_data = simulation(N, P, aP, a, b, c, d, eff.size = 2)
x = sim_data[[1]][1:round(N*r/(1+r)), , drop=F]
y = sim_data[[2]][1:round(N*r/(1+r)), , drop=F]
x.test = sim_data[[1]][(round(N*r/(1+r)) + 1):N, , drop=F]
y.test = sim_data[[2]][(round(N*r/(1+r)) + 1):N, , drop=F]

## ---------------------------------------------- ##
n_sim = 100
res_1se <- matrix(nrow = 0, ncol = 8)
res_min <- matrix(nrow = 0, ncol = 8)

for (i in seq(n_sim)){
  sim_data = simulation(N, P, aP, a, b, c, d, eff.size = 2, eff="fixed")
  x = sim_data[[1]][1:round(N*r/(1+r)), , drop=F]
  y = sim_data[[2]][1:round(N*r/(1+r)), , drop=F]
  x.test = sim_data[[1]][(round(N*r/(1+r)) + 1):N, , drop=F]
  y.test = sim_data[[2]][(round(N*r/(1+r)) + 1):N, , drop=F]

  res <- comp(x, y, x.test, y.test, Lmatrix)
  res_1se <- rbind(res_1se, res$'1se')
  res_min <- rbind(res_min, res$min)
}

meank <- function(res, type = NULL, file = NULL){
  tab_mean <- matrix(nrow=0, ncol=4)
  for(i in seq(max(res[,"k"]))) {
    tab_temp <- res[res[,"k"] == i, 1:4, drop=F]
    tab_mean <- rbind(tab_mean, colMeans(tab_temp))
  }
  tab_mean <- data.frame(tab_mean[,1:3]/tab_mean[,4])
  if(!is.null(file)) {
    jpeg(filename = file)
    matplot(tab_mean, type = c("b"),pch=1,col = c("black", "green", "red"), main = paste("SSR on number of features",type), xlab = "model size", ylab = "SSR ratio")
    legend("topright", legend = names(tab_mean),col = c("black", "green", "red"), pch=1)
    abline(h=1)
    dev.off()
  } else {
    matplot(tab_mean, type = c("b"),pch=1,col = c("black", "green", "red"), main = paste("SSR on number of features",type), xlab = "model size", ylab = "SSR ratio")
    legend("topleft", legend = names(tab_mean),col = c("black", "green", "red"), pch=1)
    abline(h=1)
  }
}

# jpeg(filename = "~/Dropbox/My R Code/VBML/comp_min.jpeg")
# boxplot(res_min[,1:4], main = "Sums of Squared Residules (min)")
# dev.off()
# jpeg(filename = "~/Dropbox/My R Code/VBML/comp_1se.jpeg")
# boxplot(res_1se[,1:4], main = "Sums of Squared Residules (1se)")
# dev.off()
# meank(res_min, type="(min)", file="~/Dropbox/My R Code/VBML/RSSk_min.jpeg")
# meank(res_1se, type="(1se)", file="~/Dropbox/My R Code/VBML/RSSk_1se.jpeg")


png(filename = "~/Dropbox/My R Code/VBML/comp_fixed.png",width = 680, height = 680)
par(mfrow=c(2,2))
boxplot(res_min[,1:5], main = "SSR (min)")
boxplot(res_1se[,1:5], main = "SSR (1se)")
meank(res_min, type="(min)")
meank(res_1se, type="(1se)")
dev.off()
