#Sim data from N(0, Omega), Omega=Lambda Lambda^T +Sigma
sim_data <- function(n, K, P, no.duke, discrete_col= NULL, alpha0, alpha_1, alpha_2){
  phi= 1 #controls error for t-dist lambdas
  
  #Set empty matrices for loops
  lambda_d <- array(0, dim=c(P, K))
  lambda_n <- array(0, dim=c(P, K))
  X <- array(0, dim=c(n, P))
  
  pi_0 <- LaplacesDemon::rdirichlet(1, rep(alpha0/K, K))
  #pi_d <- c(rep(.15, 2), runif(K-5, 0, .01), .24,  rep(.2, 2))
  pi_d <- rgamma(length(pi_0), as.vector(alpha_1*pi_0), 1)#rdirichlet(1, as.vector(alpha_1*pi_0))
  pi_n <- rgamma(length(pi_0), as.vector(alpha_2*pi_0), 1)#rdirichlet(1, as.vector(alpha_2*pi_0))
  #pi_n <- c(runif(K/2, 0, .005), .24, rep(.15, 2),runif(K/2-5, 0, .03),  rep(.19, 2))
  for(p in 2:P){
    lambda_d[p,] <- rmvnorm(1, rep(0, K), diag(as.vector(pi_d))*phi)
    lambda_n[p,] <- rmvnorm(1, rep(0, K), diag(as.vector(pi_n))*phi)
  }
  #lambda_d[1,sample(1:K, 2)] <- c(-1,1)
  #lambda_n[1,sample(1:K, 2)] <- c(-1,1)
  sigma2_d  <- rep(1, P)
  sigma2_n  <- rep(1, P)
  Omega_d <- lambda_d%*%t(lambda_d) + diag(sigma2_d)
  Omega_n <- lambda_n%*%t(lambda_n) + diag(sigma2_n)
  beta_d <- solve(Omega_d[2:P,(2:P)]) %*% Omega_d[(1), 2:P]
  beta_n <- solve(Omega_n[2:P,(2:P)]) %*% Omega_n[(1), 2:P]
  
  Duke <- sample(1:n, no.duke)
  #add in some bias
  bias <- rnorm(P, -1, 1)
  X[Duke,] <- rmvnorm(no.duke, bias, Omega_d)
  X[-Duke,] <- rmvnorm(n-no.duke, bias, Omega_n)
  
  if(length(discrete_col)>0){
    prob_x <- X[,1]
    for(j in discrete_col){
      prob_x[Duke] <- pnorm(X[Duke,j])
      prob_x[-Duke] <-  pnorm(X[-Duke,j])
      
      X[,j] <- rbinom(n, 1, prob_x)
    }
  }
  return(list(X=X, Duke=Duke,lambda_d=cbind(lambda_d, bias), lambda_n=cbind(lambda_n, bias), 
              sigma2=sigma2_d,Omega_d=Omega_d, Omega_n=Omega_n, 
              pi_0=pi_0, pi_d=pi_d, pi_n=pi_n, beta_d=beta_d, beta_n=beta_n))
}
