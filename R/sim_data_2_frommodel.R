#Sim data from X= F Lambda + Sigma
sim_data2 <- function(n, K, P, no.duke, lam=2, discrete_col= NULL, alpha0, alpha_1, alpha_2){
  
  #Set empty matrices for loops
  lambda_d <- array(0, dim=c(P, K))
  lambda_n <- array(0, dim=c(P, K))
  X <- array(0, dim=c(n, P))
  taus <- array(rexp(n*K, lam^2), dim=c(n,K))
  
  f <- array(0, dim=c(n, K))
  for(i in 1:n){
    f[i,] <- rmvnorm(1, rep(0, K), diag(taus[i,]^2))
  }
  f <- cbind(f, 1)
  pi_0 <- LaplacesDemon::rdirichlet(1, rep(alpha0/K, K))
  pi_d <- rgamma(length(pi_0), as.vector(alpha_1*pi_0), 1)#rdirichlet(1, as.vector(alpha_1*pi_0))
  pi_n <- rgamma(length(pi_0), as.vector(alpha_2*pi_0), 1)#rdirichlet(1, as.vector(alpha_2*pi_0))
  
  phi <- matrix(rgamma(P*K, 2, 2), nrow=P, ncol=K) #controls error for t-dist lambdas (tau/2)
  
  
  for(p in 1:P){
    lambda_d[p,] <- rmvnorm(1, rep(0, K), diag(as.vector(pi_d)/phi[p,]))
    lambda_n[p,] <- rmvnorm(1, rep(0, K), diag(as.vector(pi_n)/phi[p,]))
  }
  
  #lambda_d[1,sample(1:K, )] <- c(-1,1)
  #lambda_n[1,sample(1:K, 2)] <- c(-1,1)
  sigma2_d  <- rep(1, P)
  sigma2_n  <- rep(1, P)
  
  Duke <- sample(1:n, no.duke)
  #add in some bias
  bias <- rnorm(P, -3, 3)
  lambda_d <- cbind(lambda_d, bias)
  lambda_n <- cbind(lambda_n, bias)
  for(i in 1:n){
    if(i %in% Duke){
      X[i,] <- rmvnorm(1, f[i,]%*%t(lambda_d), diag(sigma2_d))
    }else{
      X[i,] <- rmvnorm(1, f[i,]%*%t(lambda_n), diag(sigma2_n))
      
    }
  }
  
  if(length(discrete_col)>0){
    prob_x <- X[,1]
    for(j in discrete_col){
      prob_x[Duke] <- pnorm(X[Duke,j])
      prob_x[-Duke] <-  pnorm(X[-Duke,j])
      
      X[,j] <- rbinom(n, 1, prob_x)
    }
  }
  return(list(X=X, Duke=Duke,lambda_d=(lambda_d), lambda_n=(lambda_n), 
              sigma2=sigma2_d,f=f, phi=phi,
              pi_0=pi_0, pi_d=pi_d, pi_n=pi_n))
}
