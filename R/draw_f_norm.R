
#Draw factors from posterior distribution under normal prior
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
