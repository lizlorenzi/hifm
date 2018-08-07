


library(mvtnorm)
library(GIGrvg)
library(Rtsne)
library(glmnet)
library(LaplacesDemon)
library(truncnorm)
library(parallel)
#library(precrec)
library(clue)


#Take posterior samples and calculate the covariance/correlation matrix
learn_covariance <- function(mod, thin){
  post_mean <- apply(mod$lambdas[thin,,,], c(2,3,4), mean)
  sigma_mean <- apply(mod$sigma2_j[thin,,], c(2,3), mean)
  K <- ncol(post_mean)-1
  P <- nrow(post_mean)
  P_Y = mod$params$P_Y
  
  omegas <- array(NA, dim=c(P, P, dim(post_mean)[3]))
  betas <- array(NA, dim=c(P-P_Y, dim(post_mean)[3]))
  for(j in 1:dim(post_mean)[3]){
    omegas[,,j] <- post_mean[,1:K,j] %*% t(post_mean[,1:K,j]) + diag(sigma_mean[,j])
    betas[,j] <- solve(omegas[(P_Y+1):P,(P_Y+1):P,j]) %*% omegas[(P_Y+1):P, (1:P_Y),j]
  }
  return(list(omegas=omegas, betas=betas))
}


#Individual parameter sampler functions
zscores <- function(y,ties.method="average"){
  z <- qnorm(rank(y,na.last="keep",ties.method=ties.method)/(sum(!is.na(y))+1))
  names(z) <- names(y)
  m <- dim(y)
  if(length(m)==2){z<-matrix(z,nrow=m[1],ncol=m[2]); dimnames(z)<-dimnames(y)}
  if(length(m)>=3){z<-array(z,dim=m); dimnames(z)<-dimnames(y)}
  return(z)
}



mvrnormR <- function(n, mu, sigma) {
  ncols <- ncol(sigma)
  mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
  mu + crossprod(t(matrix(rnorm(n * ncols), ncol = ncols)), chol(sigma))
}

draw_f_norm <- function(cut, params){
  gr = which(params$trcuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  pgroups <- params$tr_groups[gr]
  #pDuke = which(gr %in% params$oDuke)
  pxs <- params$xs[gr,] #needed for draw_xstar
  
  new_fi <- array(0, dim=c(pn, params$K))
  for(j in 1:params$J){
    pop_j <- which(pgroups==j)
    C_f <- solve(diag(params$K)+t(params$lambda_j[,1:params$K,j])%*%diag(1/params$sigma2_j[,j])%*%params$lambda_j[,1:params$K,j])
    new_fi[pop_j,] <- t(apply(pxs[pop_j,], 1, function(x) rmvnorm(1, (x%*%diag(1/params$sigma2_j[,j]) %*%(params$lambda_j[,1:params$K,j]))%*%C_f, C_f)))
    
  }
  
  #add in intercept
  fs <- cbind(new_fi, rep(1, nrow(new_fi)))
  return(fs)
}


#Draw f from laplace columnwise 
draw_fk_lap <- function( params){
  f_k <- array(0, dim=c(params$trN, params$K))
  
  for(j in 1:params$J){
    #Learn E (residuals)
    #grab intercept
    b <- tcrossprod(params$fs[which(params$tr_groups==j),params$K+1], params$lambda_j[,params$K+1,j])
    E <- params$xs[which(params$tr_groups==j),] - tcrossprod(params$fs[which(params$tr_groups==j),], params$lambda_j[,,j]) #remove info from factors 
    
    for(k in sample(1:params$K)){ #permute the draws
      #Now adjust E to account for the info this factor adds
      E <- E + tcrossprod(params$fs[which(params$tr_groups==j),k], params$lambda_j[,k,j])
      
      #Mean = XSigma^-1 Lambda_k (1+lambda_k^T Sigma^-1 Lambda_k)^-1
      mu_f <- tcrossprod(tcrossprod(E, diag(1/params$sigma2_j[,j])), t(params$lambda_j[,k,j]))
      sigma_f <- 1/((params$tau_inv[which(params$tr_groups==j),k])+c(tcrossprod(crossprod(params$lambda_j[,k,j], diag(1/params$sigma2_j[,j])), t(params$lambda_j[,k,j]))))
      f_k[which(params$tr_groups==j),k] <- rnorm(sum(params$tr_groups==j), mu_f*sigma_f, sqrt(sigma_f)+1e-8)
      E <- E - tcrossprod(f_k[which(params$tr_groups==j),k], params$lambda[,k,j])
    }
    if(any(is.na(f_k[,k]))){browser()}
    
  }
  
  fs <- cbind(f_k, rep(1, nrow(f_k)))
  return(fs)
}

inverse_gaussian <- function(mu, lambda, nsample){
  chisq1 = rnorm(nsample)^2;
  
  mu_ = mu*chisq1;
  
  r = mu + ( 0.5/lambda )*mu*( mu_ - sqrt( 4*lambda*mu_ + mu^2*chisq1^2 ) )
  
  u = runif( nsample ) 
  grab_u = u>= mu/( mu + r )
  r[grab_u] = mu[grab_u]^2/r[grab_u]
  
  return(r)
}

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


draw_xstar <- function(cut, params){
  gr = which(params$trcuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  #pDuke = which(gr %in% params$oDuke)
  pgroups <- params$tr_groups[gr]
  poX <- params$oX[gr,] #needed for draw_xstar
  pfs = params$fs[gr,]
  
  x_star <- matrix(NA, nrow=pn, ncol=params$P)  
  for(p in 1:params$P){
    if(!all(poX[,p] %in% c(0,1))){
      x_star[,p] <- poX[,p]
    }else{
      for(j in 1:params$J){
        gr_j <- which(pgroups==j)
        if(length(gr_j)>0){
          low.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,0,-Inf))
          high.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,Inf, 0))
          
          x_star[gr_j,p] <- as.vector(rtruncnorm(length(gr_j), mean=pfs[gr_j,]%*%params$lambda_j[p,,j], sd=rep(1, length(low.x)), a=low.x, b=high.x))
          
        }
      }
    }
  }
  return(x_star)
}

draw_f_test_lap <- function(cut,params,  xtest, lambda){
  gr = which(params$tecuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  pgroups <- params$te_groups[gr]
  #pDuke = which(gr %in% params$oDuke)
  pxs <- xtest[gr,] #needed for draw_xstar
  teJ <- length(unique(pgroups))
  ptaus <- params$tau_inv_test[gr,]
  
  new_fi <- array(0, dim=c(length(gr), params$K))
  # for(i in 1:params$n){ #need to run for both train and test
  #   if(i %in% params$oDuke){
  #     lambda_j <- params$lambda_d
  #     Sigma_j <- params$Sigma_d
  #   }else{
  #     lambda_j <- params$lambda_n
  #     Sigma_j <- params$Sigma_n
  #   }
  #   C_f <- solve(diag(params$K)+t(lambda_j[,1:params$K])%*%Sigma_j%*%lambda_j[,1:params$K])
  #   new_fi[i,] <- rmvnorm(1, (params$allX[i,]%*%Sigma_j %*%(lambda_j[,1:params$K]))%*%C_f, C_f)
  #only want f's to be generated from 
  # 
  
  for(j in 1:params$J){
    pop_j <- which(pgroups==j)
    for(i in pop_j){
      C_f <- solve(diag(ptaus[i,])+t(lambda[,1:params$K,j])%*%diag(1/params$sigma2_j[(params$P_Y+1):params$P,j])%*%lambda[,1:params$K,j])
      new_fi[i,] <- mvrnormR(1, (pxs[i,]%*%diag(1/params$sigma2_j[(params$P_Y+1):params$P,j]) %*%(lambda[,1:params$K,j]))%*%C_f, C_f)
    }
    
  }
  
  
  #add in intercept
  fs <- cbind(new_fi, rep(1, nrow(new_fi)))
  return(fs)
}


#Draw f from laplace columnwise 
draw_fk_test_lap <- function( params, xtest, lambda, old_fs){
  f_k <- array(0, dim=c(length(params$te_groups), params$K))
  for(j in 1:params$J){
    #Learn E (residuals)
    #grab intercept
    b <- tcrossprod(old_fs[which(params$te_groups==j),params$K+1], lambda[,params$K+1,j])
    E <- xtest[which(params$te_groups==j),] - tcrossprod(old_fs[which(params$te_groups==j),], lambda[,,j])  #remove info from factors 
    
    
    for(k in sample(1:params$K)){ #permute the draws
      #Now adjust E to account for the info this factor adds
      E <- E + tcrossprod(old_fs[which(params$te_groups==j),k], lambda[,k,j])
      
      #Mean = XSigma^-1 Lambda_k (1+lambda_k^T Sigma^-1 Lambda_k)^-1
      mu_f <- tcrossprod(tcrossprod(E, diag(1/params$sigma2_j[(params$P_Y+1):params$P,j])), t(lambda[,k,j]))
      sigma_f <- 1/((params$tau_inv_test[which(params$te_groups==j),k])+
                      tcrossprod(crossprod(lambda[,k,j], diag(1/params$sigma2_j[(params$P_Y+1):params$P,j])), t(lambda[,k,j])))
      f_k[which(params$te_groups==j),k] <- rnorm(sum(params$te_groups==j), mu_f*sigma_f, sqrt(sigma_f))
      E <- E - tcrossprod(f_k[which(params$te_groups==j),k], lambda[,k,j])
    }
    if(any(is.na(f_k[,k]))){browser()}
  }
  
  fs <- cbind(f_k, rep(1, nrow(f_k)))
  return(fs)
}

draw_f_test_norm <- function(cut,params,  xtest, lambda){
  gr = which(params$tecuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  pgroups <- params$te_groups[gr]
  #pDuke = which(gr %in% params$oDuke)
  pxs <- xtest[gr,] #needed for draw_xstar
  teJ <- length(unique(pgroups))
  
  new_fi <- array(0, dim=c(length(gr), params$K))
  # for(i in 1:params$n){ #need to run for both train and test
  #   if(i %in% params$oDuke){
  #     lambda_j <- params$lambda_d
  #     Sigma_j <- params$Sigma_d
  #   }else{
  #     lambda_j <- params$lambda_n
  #     Sigma_j <- params$Sigma_n
  #   }
  #   C_f <- solve(diag(params$K)+t(lambda_j[,1:params$K])%*%Sigma_j%*%lambda_j[,1:params$K])
  #   new_fi[i,] <- rmvnorm(1, (params$allX[i,]%*%Sigma_j %*%(lambda_j[,1:params$K]))%*%C_f, C_f)
  #only want f's to be generated from 
  # 
  
  for(j in 1:teJ){
    pop_j <- which(pgroups==j)
    C_f <- solve(diag(params$K)+t(lambda[,1:params$K,j])%*%diag(1/params$sigma2_j[(params$P_Y+1):params$P,j])%*%lambda[,1:params$K,j])
    new_fi[pop_j,] <- t(apply(pxs[pop_j,], 1, function(x) rmvnorm(1, (x%*%diag(1/params$sigma2_j[(params$P_Y+1):params$P,j]) %*%(lambda[,1:params$K,j]))%*%C_f, C_f)))
    
  }
  
  #add in intercept
  fs <- cbind(new_fi, rep(1, nrow(new_fi)))
  return(fs)
}



draw_xstar_test <- function(cut, params, lambda, ftest, xtest){
  gr = which(params$tecuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  #pDuke = which(gr %in% params$oDuke)
  pgroups <- params$te_groups[gr]
  poX <- xtest[gr,] #needed for draw_xstar
  pfs = ftest[gr,]
  teJ <- unique(pgroups)
  
  x_star <- matrix(NA, nrow=pn, ncol=params$P_X)  
  for(p in 1:(params$P_X)){
    if(!any(poX[,p] %in% c(0,1))){
      x_star[,p] <- poX[,p]
    }else{
      for(j in teJ){
        gr_j <- which(pgroups==j)
        if(length(gr_j)>0){
          
          low.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,0,-Inf))
          high.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,Inf, 0))
          x_star[gr_j,p] <- as.vector(rtruncnorm(length(gr_j), mean=pfs[gr_j,]%*%lambda[p,,j], sd=rep(1, length(low.x)), a=low.x, b=high.x))
          
        }
      }
    }
  }
  return(x_star)
}


# draw_sigma2 <- function(params){
#   new_sigma2 <- array(0, dim=c(params$P, params$J))
#   for(p in 1:params$P){
#     for(j in 1:params$J){
#       diff_j <- t(params$xs[which(params$tr_groups==j), p]-params$fs[which(params$tr_groups==j), ]%*%params$lambda_j[p,,j])%*%
#         (params$xs[which(params$tr_groups==j), p]-params$fs[which(params$tr_groups==j), ]%*%params$lambda_j[p,,j])
#       new_sigma2[p,j] <- 1/rgamma(1,params$v/2 + (sum(params$tr_groups==j))/2, 1/2*(params$v + diff_j))
#     }
#   }
#   return(new_sigma2)
# }

draw_sigma2 <- function(params){
  new_sigma2 <- array(0, dim=c(params$P, params$J))
  for(j in 1:params$J){
    diff_j <- crossprod(params$xs[which(params$tr_groups==j), ]-tcrossprod(params$fs[which(params$tr_groups==j), ],(params$lambda_j[,,j])))
    #(params$xs[which(params$tr_groups==j), params$cont_col]-params$fs[which(params$tr_groups==j), ]%*%t(params$lambda_j[params$cont_col,,j])))
    #for(p in 1:params$P){
    new_sigma2[,j] <- 1/rgamma(params$P,params$a + (sum(params$tr_groups==j))/2, 1/2*(params$b*2 + diag(diff_j)))
  }
  return(new_sigma2)
}


# draw_sigma2 <- function(params){
#   new_sigma2 <- array(0, dim=c(length(params$cont_col), params$J))
#   for(j in 1:params$J){
#     diff_j <- crossprod(params$xs[which(params$tr_groups==j), params$cont_col]-tcrossprod(params$fs[which(params$tr_groups==j), ],(params$lambda_j[params$cont_col,,j])))
#     #(params$xs[which(params$tr_groups==j), params$cont_col]-params$fs[which(params$tr_groups==j), ]%*%t(params$lambda_j[params$cont_col,,j])))
#     #for(p in 1:params$P){
#     new_sigma2[,j] <- 1/rgamma(length(params$cont_col),params$a + (sum(params$tr_groups==j))/2, 1/2*(params$b*2 + diag(diff_j)))
#   }
#   return(new_sigma2)
# }


#Draw lambda by the columns instead of rows (helps when P>>K)
draw_lambdaK <- function(params){
  #new_lambda <- array(0, dim=c(params$P, params$K+1, params$J))
  upd_lambda <- params$lambda_j #update this when we adjust using E
  for(j in 1:params$J){
    #keep the bias term 
    b <- tcrossprod(params$fs[which(params$tr_groups==j),params$K+1], upd_lambda[,params$K+1,j])
    E <- params$xs[which(params$tr_groups==j),] - tcrossprod(params$fs[which(params$tr_groups==j),], upd_lambda[,,j])  #remove info from factors 
    
    w_j <- c(params$w_j[,j],1) #to handle the intercept column
    for(k in sample(1:(params$K))){
      E <- E + tcrossprod(params$fs[which(params$tr_groups==j),k], upd_lambda[,k,j])
      
      D_lamb_inv <- diag(params$phi[,k]/w_j[k])
      sig_lamb <- solve(diag(crossprod(params$fs[which(params$tr_groups==j),k])/params$sigma2[,j]) + D_lamb_inv) #sigmas cancel out in mu (matching ricardos notation)
      mu_lamb <- crossprod(sig_lamb,crossprod(diag(params$sigma2[,j]), crossprod(E, (params$fs[which(params$tr_groups==j),k]))))
      
      upd_lambda[,k,j] <- mvrnormR(1, mu_lamb, sig_lamb)
      E <- E - tcrossprod(params$fs[which(params$tr_groups==j),k], upd_lambda[,k,j])
      
    }
    #Update the bias term
    D_lamb_inv <- diag(params$phi[,params$K+1])
    sig_lamb <- solve(diag(crossprod(params$fs[which(params$tr_groups==j),k])/params$sigma2[,j])+ D_lamb_inv) #sigmas cancel out in mu
    # <- tcrossprod(diag(1/params$sigma2[,j]), t(B))
    mu_lamb <-crossprod(sig_lamb,crossprod(diag(params$sigma2[,j]), crossprod(E, (params$fs[which(params$tr_groups==j),params$K+1]))))
    
    upd_lambda[,params$K+1,j] <- mvrnormR(1, mu_lamb, sig_lamb)
  }
  return(upd_lambda)
}

#without doing columnwise
draw_lambda <- function(params){
  new_lambda <- array(0, dim=c(params$P, params$K+1, params$J))
  for(j in 1:params$J){
    ftf <- crossprod((params$fs[which(params$tr_groups==j),]), params$fs[which(params$tr_groups==j),])
    ftx <- crossprod((params$fs[which(params$tr_groups==j), ]), params$xs[which(params$tr_groups==j),])
    for(p in 1:params$P){
      pW_inv <- diag(params$phi[p,]/c(params$w_j[,j],1) + .00001) #combine these here
      C_ld <- chol2inv(chol(pW_inv + 1/params$sigma2_j[p,j]* ftf))
      new_lambda[p,,j] <- mvrnormR(1, crossprod(t((1/params$sigma2_j[p,j]*t(ftx[,p]))), C_ld), C_ld)
    }
  }
  return(new_lambda)
}

draw_w <- function(params){
  new_w <- array(0, dim=c(params$K, params$J))
  for(j in 1:params$J){
    for(k in 1:params$K){
      Phi_k <- diag(params$phi[,k])
      new_w[k,j] <- rgig(1, lambda=params$alpha_j[j]*params$pi_0[k] - params$P/2, chi=t(params$lambda_j[,k,j])%*%Phi_k%*%params$lambda_j[,k,j], psi=2)
    }
  }
  return(new_w)
}

draw_phi <- function(params){
  phi <- array(0, dim=c(params$P, params$K+1))
  for(k in 1:(params$K+1)){
    #sum across populations for b parameter - results in P vector
    vect_B <- params$tau/2 
    for(j in 1:params$J){
      vect_B <- vect_B+ (params$lambda_j[,k,j]^2/(2*c(params$w_j[,j],1)[k]))
      
    }
    phi[,k] <- rgamma(params$P, rep(params$tau/2 + params$J/2, params$P), 
                      vect_B)
  }
  return(phi)
}




#Does draws using norm proposal - multivariate normal
draw_pi0_norm <- function(params, C=2){
  new_pi0 <- c()
  pi0 <- c()
  acc_rat <- c()
  acc <- 0
  eps <- 1e-8
  
  #location <- log(params$pi_0^2 / sqrt(C^2 + params$pi_0^2)) #need to adapt parameters so that location is at the locaiton we need
  #scale <-  sqrt(log(1 + (C^2 / params$pi_0^2)))
  #new_pi0 <- exp(rnorm(params$K,log(params$pi_0), C))
  diag(params$pi0_Sigma) = diag(params$pi0_Sigma) + 1e-8
  new_pi0 <- exp(rmvnorm(1,log(params$pi_0), C*params$pi0_Sigma))
  new_pi0 <- (new_pi0+eps)/sum(new_pi0)
  
  #draw proposal
  acc_rat_num <- ddirichlet(x=new_pi0, rep(params$alpha0/params$K, params$K), log=TRUE)+
    sum(dmvnorm(x=log(params$pi_0), log(new_pi0),  C*params$pi0_Sigma, log=TRUE))
  acc_rat_den <- ddirichlet(x=(params$pi_0), rep(params$alpha0/params$K, params$K), log=TRUE) +
    sum(dmvnorm(x=log(new_pi0), log(params$pi_0), C*params$pi0_Sigma, log=TRUE))
  
  #add on the factor specific values
  for(j in 1:params$J){
    acc_rat_num <- acc_rat_num + sum(dgamma(x=params$w_j[,j], shape=(params$alpha_j[j]*new_pi0),rate=rep(1, params$K), log=TRUE))
    acc_rat_den <- acc_rat_den + sum(dgamma(x=params$w_j[,j], shape=params$alpha_j[j]*params$pi_0, rate=rep(1, params$K), log=TRUE))
  }
  
  acc_rat <- acc_rat_num - acc_rat_den
  #print(exp(acc_rat))
  if(exp(acc_rat)>runif(1)){
    pi0 <- new_pi0
    acc=acc+1
    #if(new_pi0[k]<0){ #also reject if tries to go negative
    # pi0[k] <- params$pi_0[k]
    #  acc[k]=acc[k]-1 #adjust 
    #}
  }else{
    pi0 <- params$pi_0
  }
  return(list(pi0=pi0, acc=acc))
}

#Does draws using gamma proposal 
draw_pi0_gam <- function(params, C=1){
  new_pi0 <- c()
  pi0 <- c()
  acc_rat <- c()
  acc <- 0
  eps <- 1e-8
  
  new_theta0 <- rgamma(length(params$theta0), C*params$theta0, rate=C)
  new_pi0 <- (new_theta0+eps)/sum(new_theta0)
  
  #draw proposal
  acc_rat_num <- ddirichlet(x=new_pi0, rep(params$alpha0/params$K, params$K), log=TRUE)+
    ddirichlet(x=params$pi_0, C*new_pi0, log=TRUE)
  acc_rat_den <- ddirichlet(x=(params$pi_0), rep(params$alpha0/params$K, params$K), log=TRUE) +
    ddirichlet(x=new_pi0, C*params$pi_0, log=TRUE)
  
  #add on the factor specific values
  for(j in 1:params$J){
    acc_rat_num <- acc_rat_num + sum(dgamma(x=params$w_j[,j], shape=(params$alpha_j[j]*new_pi0),rate=rep(1, params$K), log=TRUE))
    acc_rat_den <- acc_rat_den + sum(dgamma(x=params$w_j[,j], shape=params$alpha_j[j]*params$pi_0, rate=rep(1, params$K), log=TRUE))
  }
  
  acc_rat <- acc_rat_num - acc_rat_den
  #print(exp(acc_rat))
  if(exp(acc_rat)>runif(1)){
    pi0 <- new_pi0
    theta0 <- new_theta0
    acc=acc+1
    #if(new_pi0[k]<0){ #also reject if tries to go negative
    # pi0[k] <- params$pi_0[k]
    #  acc[k]=acc[k]-1 #adjust 
    #}
  }else{
    pi0 <- params$pi_0
    theta0 <- params$theta0
  }
  return(list(pi0=pi0, acc=acc, theta0=theta0))
}



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


#Function to allign two loadings matrices - 
#Args: provide two loadings matrices (lambda1 (learned from model)  and lambda2 (truth) - the first gets reordered),
#and w_j to use when the number of factors are not the same (use to remove the zeroed columns)
reorder_loadings <- function(lambda1, lambda2, w_j=NULL){
  if(ncol(lambda1) != ncol(lambda2)){ #add on empty columns in the smaller one
    lambda2 <- cbind(lambda2, matrix(0, nrow=nrow(lambda2), ncol=(ncol(lambda1)-ncol(lambda2))))
  }
  #ordered_ws <- order(w_j, decreasing=TRUE) #find the largest weights
  #keep_cols <- ordered_ws[1:(ncol(lambda2)-1)] #keep the columns that are nonzero
  #lambda1 <- lambda1[,c(keep_cols, ncol(lambda1))] #keep nonzero columns and the intercept columns
  #}else{
  #  keep_cols=c(1:c(ncol(lambda2)-1))
  #}
  #Find distance
  dist <- crossprod(lambda2, lambda1)^2
  dist <- max(dist) - dist
  #Runs hungarian algorithm
  orders <- solve_LSAP(dist)
  print(orders)
  return(list(orders=orders, lambda1=lambda1[,orders], lambda2=lambda2, w_j = w_j[orders]))
}


#optimal lambda
opt_lambda <- function(taus){
  p <- ncol(taus)
  return(sqrt(2*p/(rowSums(1/taus))))
}


#Model sampler

#Runs Sparse TL Latent Factor model - using HDP prior on loadings matrix
#Arguments: X - predictors (unscaled), Y - outcomes, K - number of factors,
#groups - vector indicating which groups people are in, n.sim - number of iterations in sampler,
#alpha0 - pi^0 concentration parameter (for top level of HDP), alpha_j- concentration parameter should be of length J,
#v - degree of freedom for Sigma, #tau - degree of freedom for phi,
#pi0 - initial values for pi0, C - tuning parameter for pi^0 MH
#test - test indices, train-indices for training set
#fact_prior = "laplace" or "normal"
#fact_prior = "laplace" or "normal"
par_sparse_factor_model <- function(X, Y=NULL, K, groups, n.sim, alpha0, test, train, lam, fact_prior="laplace",
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


par_sparse_factor_model_wcoef <- function(X, Y=NULL, K, groups, n.sim, alpha0, test, train, lam, fact_prior="laplace",
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
  beta_j <- array(0, dim=c(n.sim, P-P_Y+1, J))
  
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
    if(it%%100==0){
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
      omega_j <- params$lambda_j[,1:K,j] %*% t(params$lambda_j[,1:K,j]) + diag(sigma2_j[it,,j])
      beta_j[it,,j] <- c(sum(params$lambda_j[,K+1,j]), solve(omega_j[(P_Y+1):P,(P_Y+1):P]) %*% omega_j[(P_Y+1):P, (1:P_Y)])
      
    }
    
    
    #if(it %% 50  == 0){
    #  print(auc(Ytest[,1], pred_resp_test[it,,1]))
    #print(head(c(Ytest[,1], pred_resp_test[it,,1])))
    #}
    
    #pred_resp_test[it,,] <- (params$fs[test,] %*% t(params$lambda_j[,,1]))[,1:P_Y]
    
    
  }
  print(it)
  #print(acceptance/n.sim)
  return(list(params=params, ftest=ftest, xtest=xtest, w_j=w_j, sigma2_j=sigma2_j, lambdas = lambdas, tau_inv_test=tau_inv_test,
              pi0=pi_0, phi=phi, acceptance=acceptance/n.sim, pred_resp = pred_resp_test, test=test, Y=Y, beta_j=beta_j))
  
}

