##
## CALCULATE GAMMA0 AND BETA0 FOR SIMULATION FOR UNIVARIATE META
##
## PREAMBLE ----------
rm(list=ls())

file.sources = list.files("../Rfn/")
sapply(paste0("../Rfn/", file.sources), source)

library(dplyr)
library(MASS)
library(mnormt)
library(ggplot2)
library(foreach)
library(doRNG)
library(metafor)


## true parameters set1 ----------
s = c(10, 50) ## #of population studies
set = expand.grid(
  t.theta = c(0.7), ## true theta
  t.tau = sqrt(c(0.05, 0.3, 0.7)), ## true tau  
  t.rho = c(0.8), ## for population data rho does not matter the estimates
  n.median = c(10), ## median number of total subjects,
  grp.r = c(1, 2), ##group ratio: treat:control
  pmax = 0.99,
  pmin = 0.6
) %>% arrange(t.theta,n.median)
# set$p0 = ifelse(set$n.median==50, 0.05,0.1)
# set$t.rho = ifelse(set$t.theta>0, 0.8, -0.8)
# set$ymin = ifelse(set$n.median==50, 0, 0)
# set$ymax = ifelse(set$n.median==50, 5, 10)
set$ymin = 0
set$ymax = set$n.median/2
# set.gr1 = set[set$grp.r==1,]
# set.gr2 = set[set$grp.r==2,]
# set.gr = set.gr1[1,]

## ----------
## SIMULATION 
ncores = parallel::detectCores()-1
cl = parallel::makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

set.seed(2025)
for(S in s){
for(i in 1:nrow(set)){ 
  DATA = foreach(r=1:200, .combine = rbind,.packages=c("tidyr","mnormt","dplyr","metafor"))  %dorng%  {

    set.gr = set[i,]
    pmax = set.gr$pmax
    pmin = set.gr$pmin
    
    plist1 = gen.data1(
      s=S, n.med=set.gr$n.median, gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      y_min=set.gr$ymin,y_max=set.gr$ymax,
      Pnmax = pmax, Pnmin = pmin)

    
    plist2 = gen.data2(
      s=S, n.med=set.gr$n.median, gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      p0 = set.gr$p0,
      Pnmax=pmax,Pnmin=pmin
      )
    
    
    ## population data and selective data
    pdata1 = plist1$p.dt
    pdata2 = plist2$p.dt

    ## data of lnOR
    data = lapply(list(pdata1,pdata2),
           function(data) escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data)
             )
    lpdata1 = data[[1]]
    lpdata2 = data[[2]]

    ## set parset list
    parset.nn = list(
          mu.bound = 1, 
          tau.bound = 1,
          eps = 1e-3,
          init.vals = c(0.1,0.01)
          )
    
    ## estimation without/with PB: NN, HN-GLMM, BN-GLMM on pdata and sdata
    fit.nn = lapply(
      list(lpdata1,lpdata2), 
      function(data) with(data, NN_LMM(yi, vi, 
        parset=parset.nn)))
    
    pnn1 = c(fit.nn[[1]]$mu, fit.nn[[1]]$tau, 
      cv = ifelse(is.null(fit.nn[[1]]$opt$convergence), NA, fit.nn[[1]]$opt$convergence))
    pnn2 = c(fit.nn[[2]]$mu, fit.nn[[2]]$tau, 
      cv = ifelse(is.null(fit.nn[[2]]$opt$convergence), NA, fit.nn[[2]]$opt$convergence))

    ## initial values for GLMM
    parset.glmm = list(
      mu.bound = 1,
      tau.bound = 1.5,
      eps = 1e-3,
      integ.limit = 5, 
      cub.tol = 1e-5,
      init.vals = c(0.1,0.01)
    )
    
    
    fit.hn = lapply(
      list(lpdata1,lpdata2), 
      function(data) with(data, HN_GLMM(y0, y1, n0, n1, 
        parset = parset.glmm)))
    phn1 = c(fit.hn[[1]]$mu, fit.hn[[1]]$tau,
      cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[1]]$opt$convergence))
    phn2 = c(fit.hn[[2]]$mu, fit.hn[[2]]$tau, 
      cv = ifelse(is.null(fit.hn[[2]]$opt$convergence), NA, fit.hn[[2]]$opt$convergence))

    fit.bn = lapply(
      list(lpdata1,lpdata2), 
      function(data) with(data, BN_GLMM(y0, y1, n0, n1, 
        parset = parset.glmm)))
    pbn1 = c(fit.bn[[1]]$mu, fit.bn[[1]]$tau, 
      cv = ifelse(is.null(fit.bn[[1]]$opt$convergence), NA, fit.bn[[1]]$opt$convergence))
    pbn2 = c(fit.bn[[2]]$mu, fit.bn[[2]]$tau, 
      cv = ifelse(is.null(fit.bn[[2]]$opt$convergence), NA, fit.bn[[2]]$opt$convergence))
      
   
    
    res = rbind(
      pnn1,pnn2,
      phn1,phn2,
      pbn1,pbn2)

  }
  save(DATA,file = paste0("res-pop/data-set-",i,"-S",S,".RData"))
  
}}

parallel::stopCluster(cl)

