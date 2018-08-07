#Draw weights w - from posterior draw of generalized inverse gaussian
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
