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


## true parameters ----------
s = c(10, 25, 50) ## #of population studies
set = expand.grid(
  t.theta = c(-0.7,0.7), ## true theta
  t.tau = c(0.05, 0.15, 0.6), ## true tau  
  # t.rho = c(-0.8,0.8),
  n.median = c(50, 100), ## median number of total subjects
  grp.r = c(1, 2), ##group ratio: treat:control
  pmax = 0.99,
  pmin = 0.6
) %>% arrange(t.theta,n.median)
set$p0 = ifelse(set$n.median==50, 0.05,0.1)
set$t.rho = ifelse(set$t.theta>0, 0.8, -0.8)

set.gr1 = set[set$grp.r==1,]
set.gr2 = set[set$grp.r==2,]

## ----------
## SIMULATION 
ncores = parallel::detectCores()-1
cl = parallel::makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

set.gr = set.gr1
set.seed(2024)
for(S in s[1]){
for(i in 1){ 
  DATA = foreach(r=1:1000, .combine = rbind,.packages=c("tidyr","mnormt","dplyr"),.errorhandling="remove")  %dorng%  {
  
    S = 100
    i = 1
    pmax = set.gr$pmax[i]
    pmin = set.gr$pmin[i]
    plist = gen.data3(
      s=S,
      theta=set.gr$t.theta[i],tau=set.gr$t.tau[i],
      rho=set.gr$t.rho[i],
      p0 = set.gr$p0[i],
      n.med =set.gr$n.median[i],#set.val$n0.median,
      gr=set.gr$grp.r[i],
      Pnmax=pmax,Pnmin=pmin)
    
    
    ## population data and selective data
    pdata = plist$p.dt
    sdata = plist$s.dt
    
    ## data of lnOR
    data = lapply(list(pdata,sdata),
           function(data) escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data)
             )
    pdata = data[[1]]
    sdata = data[[2]]

    ## set parset list
    parset.nn = list(
          mu.bound = 2, 
          tau.bound = 1,
          eps = 1e-3,
          init.vals = c(set.gr$t.theta[i],set.gr$t.tau[i])
          )
    parset.glmm = list(
          mu.bound = 2,
          tau.bound = 1,
          eps = 1e-3,
          integ.limit = 10, 
          cub.tol = 1e-5,
          init.vals = c(set.gr$t.theta[i],set.gr$t.tau[i])
          )
    ## estimation without/with PB: NN, HN-GLMM, BN-GLMM on pdata and sdata
    fit.nn = lapply(
      list(pdata,sdata), 
      function(data) with(data, NN_LMM(yi, vi, 
        parset=parset.nn)))
    pnn = c(fit.nn[[1]]$mu, fit.nn[[1]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.nn[[1]]$opt$convergence), NA, fit.nn[[1]]$opt$convergence))
    snn = c(fit.nn[[2]]$mu, fit.nn[[2]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.nn[[2]]$opt$convergence), NA, fit.nn[[2]]$opt$convergence))

    fit.hn = lapply(
      list(pdata,sdata), 
      function(data) with(data, HN_GLMM(y0, y1, n0, n1, 
        parset = parset.glmm)))
    phn = c(fit.hn[[1]]$mu, fit.hn[[1]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[1]]$opt$convergence))
    shn = c(fit.hn[[2]]$mu, fit.hn[[2]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.hn[[2]]$opt$convergence), NA, fit.hn[[2]]$opt$convergence))

    fit.bn = lapply(
      list(pdata,sdata), 
      function(data) with(data, BN_GLMM(y0, y1, n0, n1, 
        parset = parset.glmm)))
    pbn = c(fit.bn[[1]]$mu, fit.bn[[1]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.bn[[1]]$opt$convergence), NA, fit.bn[[1]]$opt$convergence))
    
    sbn = c(fit.bn[[2]]$mu, fit.bn[[2]]$tau, rho=rep(NA,2), 
      cv = ifelse(is.null(fit.bn[[2]]$opt$convergence), NA, fit.bn[[2]]$opt$convergence))
      
    ## set parset list
    parset.c20 = list(
      mu.bound = 2,
      tau.bound = 1,
      estimate.rho = TRUE, 
      eps = 1e-3,
      init.vals = c(set.gr$t.theta[i],set.gr$t.tau[i])
    )
    parset.prop = list(
      mu.bound = 2,
      tau.bound = 1,
      estimate.rho = TRUE, 
      eps = 1e-3,
      integ.limit = 10, 
      cub.tol = 1e-5,
      init.vals = c(set.gr$t.theta[i],set.gr$t.tau[i])
    )
    ## adjust for PB: COPAS2000, COPAS_HNGLMM, COPAS_BNGLMM
    fit.c20 = suppressWarnings(with(sdata, COPAS2000(yi, vi, Psemax =pmin, Psemin = pmax,
      parset=parset.c20)))
    cnn = c(fit.c20$mu, fit.c20$tau, fit.c20$rho,
      cv = ifelse(is.null(fit.c20$opt$convergence), NA, fit.c20$opt$convergence))

    fit.chn = suppressWarnings(with(sdata, COPAS_HNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
      parset=parset.prop)))
    chn = c(fit.chn$mu, fit.chn$tau, fit.chn$rho,
      cv = ifelse(is.null(fit.chn$opt$convergence), NA, fit.chn$opt$convergence))
    
    fit.cbn = suppressWarnings(with(sdata, COPAS_BNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
      parset=parset.prop)))
    cbn = c(fit.cbn$mu, fit.cbn$tau, fit.cbn$rho,
      cv = ifelse(is.null(fit.cbn$opt$convergence), NA, fit.cbn$opt$convergence))

    res = rbind(
      pnn,snn,cnn,
      phn,shn,chn,
      pbn,sbn,cbn)

  }
  save(DATA,file = paste0("res/data-sen-",i,"-S",SS,".RData"))
  
}}

parallel::stopCluster(cl)

## summary
load("res/data-sen-2-S15.RData")
DATA[DATA[,8]==1,]=NA
rr = dim(DATA)[1]/9
group_indices = rep(1:9, each = rr)
rowsum(DATA, group_indices, na.rm = T)/rr
