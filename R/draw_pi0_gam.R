#Draw pi0 using gamma proposal within MH step
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
