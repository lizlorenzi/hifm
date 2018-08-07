#Draw scale parameter for loadings matrix, phi, from gamma distribution
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
