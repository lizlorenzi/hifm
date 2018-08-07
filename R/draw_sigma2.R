#Draws sigma2 from invgamma(tau/2, tau/2) prior
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
