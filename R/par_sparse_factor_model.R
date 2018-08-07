#' Runs Hierarchical infinite factor model - using HDP prior on loadings matrix
#' @param X  predictors (unscaled).
#' @param Y  outcomes.
#' @param K  number of factors.
#' @param groups  vector indicating which groups people are in (for single group, rep(1, nrow(X)))
#' @param n.sim  number of iterations in sampler.
#' @param alpha0  pi^0 concentration parameter (for top level of HDP).
#' @param alpha_j concentration parameter should be of length J.
#' @param v  degree of freedom for Sigma
#' @param tau  degree of freedom for phi.
#' @param pi0  initial values for pi0
#' @param C  tuning parameter for pi^0 MH
#' @param test  test indices
#' @param train indices for training set
#' @param fact_prior = "laplace" or "normal"
#' @return List of multiple return arguments: params - final iteration of all parameters,
#' ftest - posterior samples of factors for test set (without any information on y),
#' xtest - posterior samples of test set x with transformations, w_j - iterations of weights for each population,
#' sigma2_j - posterior iterations of sigma2 (idiosyncratic noise), lambdas - posterior of loadings matrix,
#' pred_resp - posterior predictive response of test set
#' @examples   
#' sim1 <- sim_data(500, 5, 20, 300,alpha0 = 25, alpha_1=25, alpha_2=25)
#' test <- sample(sim1$Duke, 100);
#' groups_sim <- rep(1, 500); groups_sim[-sim1$Duke]=2
#' test_norm <- hifm(sim1$X[,-c(1)], Y=sim1$X[-test,1],K= 10,
#'                                           groups=groups_sim, n.sim=1000, alpha0=15,test= test,
#'                                           train=c(1:500)[-test], alpha_j=c(20,20), a=5, b=4, tau=4,
#'                                           lam=2, J=2, C=40,fact_prior="normal")

hifm <- function(X, Y=NULL, K, groups, n.sim, alpha0, test, train, lam, fact_prior="laplace",
                                    alpha_j, a, b, tau, pi0=NULL, wj1=NULL, wj2=NULL,phi_0=NULL, C=2, J=2, Ytest=NULL){
  no_cores =no.cores= detectCores() - 1
  
  Z <- as.matrix(cbind(Y,X[train,]))
  n <- nrow(X)
  trn <- length(train)
  te_j <- length(unique(groups[test]))
  P_X <- ncol(X)
  P_Y <- ncol(Y)
  if(length(P_Y)==0){
    if(length(Y)>0){ #for one outcome - reads as vector so P_Y=NULL
      P_Y=1
    }else{
      P_Y=0
    }
  }
  P <- P_X+P_Y
  
  w_j <- array(1/K, dim=c(n.sim, K, J))
  pi_0 <- array(1/K, dim=c(n.sim, K))
  sigma2_j <- array(1, dim=c(n.sim, P, J))
  tau_inv <- array(1, dim=c(length(train), K))
  tau_inv_test <- array(1, dim=c(length(test), K))
  #cont_col <- which(apply(Z, 2, function(x) any(!x %in% c(0,1))))
  cont_col <- which(apply(Z[1:100,], 2, function(x) any(x>2)))
  
  #sigma2_j[,cont_col,] <- .1
  
  phi <- array(2, dim=c(n.sim, P, K+1))
  #don't need to store each iteration from below since non-identifiable
  lambdas <- array(0,dim=c(n.sim, P, K+1, J)) #to save 
  #Set the intercept equal to the mean of each covariate
  lambda_j <- array(rnorm(P*(K+1)*J, 0, 1), dim=c(P, K+1, J))
  for(j in 1:J){
    lambda_j[cont_col,K+1, j] = colMeans(Z[which(groups[train]==j),cont_col])
    #lambda_j[-cont_col,K+1, j] = qnorm(colMeans(Z[which(groups[train]==j),-cont_col]))
    if(any(is.infinite(lambda_j[-cont_col,K+1, j]))){
      lambda_j[-cont_col,K+1, j][is.infinite(lambda_j[-cont_col,K+1, j])] = 0
    }
  }
  
  F <- array(rnorm(n*K, 0, 1), dim=c(n, K))
  F <- cbind(F, rep(1, nrow(F)))
  pred_resp_test <- array(0, dim=c(n.sim, length(test), P_Y))
  pred_resp_train <- array(0, dim=c(n.sim, 1000, P_Y))
  #for MH step
  acceptance <- 0
  
  #Draw pi_0
  if(length(pi0)==0){
    theta0 <- rgamma(K, alpha0/K, 1)
    pi_0[1,] <- theta0/sum(theta0)  
  }else{
    if(length(pi0)<K){
      
      pi0 <- c(pi0, rep(1e-8, K-length(pi0)))
      print(pi0)
    }
    theta0 <- rgamma(K, alpha0/K, 1)
    
    pi_0[1,] <- pi0 #inits provided
    #w_j[1,,1] <- wj1
    #w_j[1,,2] <- wj2
    #phi[1,,] <- phi_0
  }
  
  #Save original data
  oX <- Z
  oXtest <- X[test,]
  
  tr_groups <- groups[train]
  te_groups <- groups[test]
  
  ftest <- array(0, dim=c(n.sim, length(test), K+1))
  ftest[1,,] <-  F[1:length(test),] #just to initialize
  #Seal all the parameters into a list
  params <- list(n=n, trN = trn, P=P,P_Y=P_Y,P_X=P_X, K=K, J=J, w_j=matrix(w_j[1,,], ncol=J), pi_0=pi_0[1,], theta0=theta0, test=test, train=train,
                 sigma2_j=matrix(sigma2_j[1,,], ncol=J), alpha_j=alpha_j, alpha0=alpha0, lambda_j=lambda_j,oXtest = oXtest, lam=lam, tau_inv=tau_inv,tau_inv_test=tau_inv_test, 
                 fs=F, F=F, a=a,b=b, phi=phi[1,,], tau=tau, groups=groups, tr_groups=tr_groups, te_groups=te_groups,allX=Z, oX=oX, xs=oX, cont_col=cont_col, te_j = length(unique(te_groups)))  
  
  for(it in 2:n.sim){
    if(it%%10==0){
      print(c(it, acceptance/it))
      
    }
    
    #First update all on only training data
    
    #Start parallelization cluster
    #cut from the training indices and all indices (need to do the xstar draws different for test data)
    params$trcuts <- cut(1:trn, no.cores, labels=(1:(no.cores)))
    #params$cuts <- cut(1:n, no_cores-1, labels=(1:(no_cores-1)))
    
    xs <- mclapply(unique(params$trcuts), draw_xstar, params, mc.cores = no.cores)
    params$xs <- do.call(rbind, xs)
    #params$xs[train,] <- params$X
    
    if(fact_prior=="laplace"){
      #fs <- mclapply(unique(params$trcuts), draw_f_lap, params, mc.cores=10)
      #params$fs <- do.call(rbind, fs)
      fs <- draw_fk_lap(params)
      params$fs <- fs
      tau_inv <- draw_tau_inv(params, params$fs)
      params$tau_inv <- tau_inv
    }else{
      fs <- draw_fk_lap(params)
      params$fs <- fs
      #fs <- mclapply(unique(params$trcuts), draw_f_norm, params, mc.cores=no_cores)
      #params$fs <- do.call(rbind, fs)
    }
    
    
    #params$F[train,] <- params$fs
    
    #Stop parallel cluster
    #stopCluster(cl)
    
    
    params$lambda_j <- draw_lambda(params)
    lambdas[it,,,] <- params$lambda_j
    #Subset lambda for just xs - used for test set
    lambda_x <- array(params$lambda_j[((P_Y+1):P),,], dim=c(P_X, K+1, J))
    
    
    
    
    sigma2_j[it,,] <- draw_sigma2(params)
    params$sigma2_j = matrix(sigma2_j[it,,], ncol=J)
    
    
    pi_0_res <- draw_pi0_gam(params, C)
    pi_0[it,] <- pi_0_res$pi0
    params$theta0 <- pi_0_res$theta0
    acceptance <- acceptance + pi_0_res$acc
    params$pi_0 <- pi_0[it,]
    phi[it,,] <- draw_phi(params)
    params$phi = phi[it,,]
    
    w_j[it,,] <- draw_w(params)
    params$w_j <- matrix(w_j[it,,], ncol=J)
    # 
    # 
    
    
    if(it%%50==0){
      plot(w_j[it,,1])
    }
    
    #Now evaluate for test data
    params$tecuts <- cut(1:length(test), no.cores-1, labels=(1:(no.cores-1)))
    xtest_l <- mclapply(unique(params$tecuts), draw_xstar_test, params, lambda_x,  ftest[it-1,,], params$oXtest, mc.cores=no.cores)
    xtest <- do.call(rbind, xtest_l)
    if(fact_prior=="laplace"){
      #ftest_l <- mclapply(unique(params$tecuts), draw_f_test_lap, params, xtest, lambda_x,mc.cores=10)
      #ftest[it,,] <- do.call(rbind, ftest_l)
      ftest[it,,] <- draw_fk_test_lap(params, xtest, lambda_x, ftest[it-1,,])
      
      tau_inv_test <- draw_tau_inv(params, ftest[it,,])
      params$tau_inv_test <- tau_inv_test
    }else{
      ftest[it,,] <- draw_fk_test_lap(params, xtest, lambda_x, ftest[it-1,,])
      
    }
    
    for(j in 1:J){
      gr_j = which(te_groups==j)
      if(length(gr_j)>0){
        pred_resp_test[it,gr_j,] <- (ftest[it,gr_j,] %*% t(matrix(params$lambda_j[1:P_Y,,j], nrow=P_Y, ncol=ncol(ftest[it,gr_j,]))))
      }
    }
    
    
    #if(it %% 50  == 0){
    #  print(auc(Ytest[,1], pred_resp_test[it,,1]))
    #print(head(c(Ytest[,1], pred_resp_test[it,,1])))
    #}
    
    #pred_resp_test[it,,] <- (params$fs[test,] %*% t(params$lambda_j[,,1]))[,1:P_Y]
    #omega_d <- params$lambda_j[,1:K,1] %*% t(params$lambda_j[,1:K,1]) + diag(sigma2_j[it,,1])
    #beta_d[it,,] <- solve(omega_d[(P_Y+1):P,(P_Y+1):P]) %*% omega_d[(P_Y+1):P, (1:P_Y)]
    
    #omega_n <- params$lambda_j[,1:K,2] %*% t(params$lambda_j[,1:K,2]) + diag(sigma2_j[it,,2])
    #beta_n[it,,] <- solve(omega_n[(P_Y+1):P,(P_Y+1):P]) %*% omega_n[(P_Y+1):P, (1:P_Y)]
  }
  print(it)
  #print(acceptance/n.sim)
  return(list(params=params, ftest=ftest, xtest=xtest, w_j=w_j, sigma2_j=sigma2_j, lambdas = lambdas, tau_inv_test=tau_inv_test,
              pi0=pi_0, phi=phi, acceptance=acceptance/n.sim, pred_resp = pred_resp_test, test=test, Y=Y))
  
}
