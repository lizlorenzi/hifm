#Draw factors (under laplace prior) for testing cohort - columnwise
#partitions lambda appropriately to avoid using y information for test set
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
