#Draw factors from posterior distribution under laplace prior - columnwise for efficiency
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
