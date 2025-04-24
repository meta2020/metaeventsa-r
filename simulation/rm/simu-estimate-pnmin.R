##
## CALCULATE GAMMA0 AND BETA0 FOR SIMULATION FOR UNIVARIATE META
##
## PREAMBLE ----------
rm(list=ls())

source("preamble.R")
## ----------
## SIMULATION 
ncores <- parallel::detectCores()-1
cl <- parallel::makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

set.seed(2024)

pnmin0 = NULL
for(ir in 1:4){ 
  
  DATA <- foreach(r=1:100, .combine = rbind,.packages=c("tidyr","mnormt","dplyr"))  %dorng%  {
    
    # ir=1
    set.val = set[ir,]
    pnmin = gen.pnmin(
      s=5000,theta=set.val$theta,tau=set.val$tau,rho=set.val$rho,
      n_min=set.val$n.min,n_max=set.val$n.max,
      y_min=set.val$y.min,y_max=set.val$n.max,
      Pnmax=set.val$pnmax,p=0.7)
    
    return(c(ir,pnmin))

  }
  
  pnmin0 = rbind(pnmin0, colMeans(DATA))
  }

write.csv(round(pnmin0, 5), "para10.csv")

parallel::stopCluster(cl)
