
#Performs probit transform for binary columns of X 
#Returns X matrix with binary columns converted
draw_xstar <- function(cut, params){
  gr = which(params$trcuts==cut)
  #Create objects in params to allow for the par.
  pn = length(gr) #new n for parallelized runs
  #pDuke = which(gr %in% params$oDuke)
  pgroups <- params$tr_groups[gr]
  poX <- params$oX[gr,] #needed for draw_xstar
  pfs = params$fs[gr,]
  
  x_star <- matrix(NA, nrow=pn, ncol=params$P)  
  for(p in 1:params$P){
    if(!all(poX[,p] %in% c(0,1))){
      x_star[,p] <- poX[,p]
    }else{
      for(j in 1:params$J){
        gr_j <- which(pgroups==j)
        if(length(gr_j)>0){
          low.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,0,-Inf))
          high.x <- sapply(poX[gr_j,p], function(x) ifelse(x==1,Inf, 0))
          
          x_star[gr_j,p] <- as.vector(rtruncnorm(length(gr_j), mean=pfs[gr_j,]%*%params$lambda_j[p,,j], sd=rep(1, length(low.x)), a=low.x, b=high.x))
          
        }
      }
    }
  }
  return(x_star)
}
