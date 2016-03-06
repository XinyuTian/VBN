pred_log = function(x, y, beta) {
  y = y / 2 + 0.5
  ypred = round(sigmoid(x %*% matrix(beta, ncol = 1)))
  return(mean(abs((ypred - y) / 2)))
}