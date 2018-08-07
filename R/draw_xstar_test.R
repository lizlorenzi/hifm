#Draw transformation of X for testing cohort - 
#partitions lambda appropriately to avoid using y information for test set
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
