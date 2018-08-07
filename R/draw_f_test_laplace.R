#Draw factors (under laplace prior) for testing cohort - 
#partitions lambda appropriately to avoid using y information for test set
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
