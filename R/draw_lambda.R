#Draw loadings matrix for each population
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
