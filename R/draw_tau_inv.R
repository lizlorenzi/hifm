
#Draw posterior for 1/tau - auxiliary parameter for laplace factor model
draw_tau_inv <- function(params, fs){
  #assume sigma2=1
  mu_p <- array(0, dim=c(nrow(fs), params$K)) #params for the inv gaussian
  lambda_p <- c()
  tau_inv_sq <- array(0, dim=c(nrow(fs), params$K))
  for(k in 1:params$K){
    mu_p[,k] <- sqrt(params$lam^2/fs[,k]^2)
    lambda_p[k] <- params$lam^2
    tau_inv_sq[,k] <- sapply(mu_p[,k], inverse_gaussian, lambda_p[k],1)
  }
  return(tau_inv_sq)
}

