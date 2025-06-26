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
  t.theta = c(-0.7), ## true theta
  t.tau = sqrt(c(0.05, 0.3, 0.7)), ## true tau  
  t.rho = c(0.8), ## for population data rho does not matter the estimates
  n.median = c(20), ## median number of total subjects,
  grp.r = c(1, 2), ##group ratio: treat:control
  pmax = 0.99,
  pmin = 0.5
) %>% arrange(t.theta,n.median)
set$p0 = 0.05
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
    
    ## create datasets
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
    
    
    ## selected data and population data
    sdata1 = plist1$s.dt
    sdata2 = plist2$s.dt
    pdata1 = plist1$p.dt
    pdata2 = plist2$p.dt
    
    ## data with n and lnOR
    data = lapply(list(pdata1,pdata2, sdata1,sdata2),
                  function(data) escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data)
    )
    lpdata1 = data[[1]]
    lpdata2 = data[[2]]
    lsdata1 = data[[3]]
    lsdata2 = data[[4]]
    
    ## population and selected models -------
    
    ## set parset list
    parset.nn = list(
      mu.bound = 1, 
      tau.bound = 1,
      eps = 1e-3,
      init.vals = c(0.1,0.01)
    )
    
    ## estimation without/with PB: NN, HN-GLMM, BN-GLMM on pdata and sdata
    fit.nn = lapply(
      list(lpdata1,lpdata2, lsdata1, lsdata2), 
      function(data) with(data, NN_LMM(yi, vi, parset=parset.nn)))
    
    pnn1 = c(fit.nn[[1]]$mu, fit.nn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[1]]$opt$convergence), NA, fit.nn[[1]]$opt$convergence))
    pnn2 = c(fit.nn[[2]]$mu, fit.nn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[2]]$opt$convergence), NA, fit.nn[[2]]$opt$convergence))
    snn1 = c(fit.nn[[3]]$mu, fit.nn[[3]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[3]]$opt$convergence), NA, fit.nn[[3]]$opt$convergence))
    snn2 = c(fit.nn[[4]]$mu, fit.nn[[4]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[4]]$opt$convergence), NA, fit.nn[[4]]$opt$convergence))
    
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
      list(lpdata1,lpdata2, lsdata1, lsdata2), 
      function(data) with(data, HN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    phn1 = c(fit.hn[[1]]$mu, fit.hn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[1]]$opt$convergence))
    phn2 = c(fit.hn[[2]]$mu, fit.hn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[2]]$opt$convergence), NA, fit.hn[[2]]$opt$convergence))
    shn1 = c(fit.hn[[3]]$mu, fit.hn[[3]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[3]]$opt$convergence))
    shn2 = c(fit.hn[[4]]$mu, fit.hn[[4]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[4]]$opt$convergence), NA, fit.hn[[4]]$opt$convergence))
    
    fit.bn = lapply(
      list(lpdata1,lpdata2, lsdata1, lsdata2), 
      function(data) with(data, BN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    pbn1 = c(fit.bn[[1]]$mu, fit.bn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[1]]$opt$convergence), NA, fit.bn[[1]]$opt$convergence))
    pbn2 = c(fit.bn[[2]]$mu, fit.bn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[2]]$opt$convergence), NA, fit.bn[[2]]$opt$convergence))
    sbn1 = c(fit.bn[[3]]$mu, fit.bn[[3]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[3]]$opt$convergence), NA, fit.bn[[3]]$opt$convergence))
    sbn2 = c(fit.bn[[4]]$mu, fit.bn[[4]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[4]]$opt$convergence), NA, fit.bn[[4]]$opt$convergence))
    
    ## adjusted models -------
    ## Copas 1990
    parset.copas = list(
    mu.bound = 1,
    tau.bound = 1,
    estimate.rho = TRUE, 
    eps = 1e-3,
    init.vals = c(0.1,0.1,0.1) ## initials for mu tau and rho
    )
    fit.C90 = lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS1999(yi, n, Pnmax = pmax, Pnmin = pmin,
                                    parset=parset.copas)))
    
    pc901 = c(fit.C90[[1]]$mu, fit.C90[[1]]$tau,fit.C90[[1]]$rho,
      cv = ifelse(is.null(fit.C90[[1]]$opt$convergence), NA, fit.C90[[1]]$opt$convergence))
    pc902 = c(fit.C90[[2]]$mu, fit.C90[[2]]$tau,fit.C90[[1]]$rho,
      cv = ifelse(is.null(fit.C90[[2]]$opt$convergence), NA, fit.C90[[2]]$opt$convergence))
    
    ## Copas 2000
    fit.C20 = lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS2000(yi, vi, Psemax = pmax, Psemin = pmin,
                                    parset=parset.copas)))
    
    pc201 = c(fit.C20[[1]]$mu, fit.C20[[1]]$tau,fit.C20[[1]]$rho,
              cv = ifelse(is.null(fit.C20[[1]]$opt$convergence), NA, fit.C20[[1]]$opt$convergence))
    pc202 = c(fit.C20[[2]]$mu, fit.C20[[2]]$tau,fit.C20[[1]]$rho,
              cv = ifelse(is.null(fit.C20[[2]]$opt$convergence), NA, fit.C20[[2]]$opt$convergence))

    
    ## New copas HNGLMM
    parset.new = list(
      mu.bound = 1,
      tau.bound = 1,
      estimate.rho = TRUE, 
      eps = 1e-3,
      integ.limit = 5, 
      cub.tol = 1e-5,
      init.vals = c(0.1,0.1,0.1)
    )
    
    
    fit.Chn = suppressWarnings(lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS_HNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                       n_min = min(sdata1$n), n_max = max(sdata1$n),
                                       parset=parset.new))))
    
    pchn1 = c(fit.Chn[[1]]$mu, fit.Chn[[1]]$tau,fit.Chn[[1]]$rho,
              cv = ifelse(is.null(fit.Chn[[1]]$opt$convergence), NA, fit.Chn[[1]]$opt$convergence))
    pchn2 = c(fit.Chn[[2]]$mu, fit.Chn[[2]]$tau,fit.Chn[[1]]$rho,
              cv = ifelse(is.null(fit.Chn[[2]]$opt$convergence), NA, fit.Chn[[2]]$opt$convergence))
    
    ## New copas HNGLMM with fixed rho
    parset.new.fr = list(
      mu.bound = 1,
      tau.bound = 1,
      estimate.rho = FALSE, 
      rho.fix = 0.8, 
      eps = 1e-3,
      integ.limit = 5, 
      cub.tol = 1e-5,
      init.vals = c(0.1,0.1)
    )
    
    fit.Chnf = suppressWarnings(lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS_HNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                       n_min = min(sdata1$n), n_max = max(sdata1$n),
                                       parset=parset.new.fr))))
    
    pchnf1 = c(fit.Chnf[[1]]$mu, fit.Chnf[[1]]$tau,fit.Chnf[[1]]$rho,
              cv = ifelse(is.null(fit.Chnf[[1]]$opt$convergence), NA, fit.Chnf[[1]]$opt$convergence))
    pchnf2 = c(fit.Chnf[[2]]$mu, fit.Chnf[[2]]$tau,fit.Chnf[[1]]$rho,
              cv = ifelse(is.null(fit.Chnf[[2]]$opt$convergence), NA, fit.Chnf[[2]]$opt$convergence))

    ## New copas BNGLMM
    fit.Cbn = suppressWarnings(lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS_BNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                       n_min = min(sdata1$n), n_max = max(sdata1$n),
                                       parset=parset.new))))
    
    pcbn1 = c(fit.Cbn[[1]]$mu, fit.Cbn[[1]]$tau,fit.Cbn[[1]]$rho,
              cv = ifelse(is.null(fit.Cbn[[1]]$opt$convergence), NA, fit.Cbn[[1]]$opt$convergence))
    pcbn2 = c(fit.Cbn[[2]]$mu, fit.Cbn[[2]]$tau,fit.Cbn[[1]]$rho,
              cv = ifelse(is.null(fit.Cbn[[2]]$opt$convergence), NA, fit.Cbn[[2]]$opt$convergence))
    
    ## New copas BNGLMM with fixed rho
    fit.Cbnf = suppressWarnings(lapply(
      list(lsdata1, lsdata2), 
      function(data) with(data, 
                          COPAS_BNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                       n_min = min(sdata1$n), n_max = max(sdata1$n),
                                       parset=parset.new.fr))))
    
    pcbnf1 = c(fit.Cbnf[[1]]$mu, fit.Cbnf[[1]]$tau,fit.Cbnf[[1]]$rho,
              cv = ifelse(is.null(fit.Cbnf[[1]]$opt$convergence), NA, fit.Cbnf[[1]]$opt$convergence))
    pcbnf2 = c(fit.Cbnf[[2]]$mu, fit.Cbnf[[2]]$tau,fit.Cbnf[[1]]$rho,
              cv = ifelse(is.null(fit.Cbnf[[2]]$opt$convergence), NA, fit.Cbnf[[2]]$opt$convergence))
   
    
    res = rbind(
      pnn1,phn1,pbn1,snn1,shn1,sbn1,pc901,pc201,pchn1,pchnf1,pcbn1,pcbnf1,
      pnn2,phn2,pbn2,snn2,shn2,sbn2,pc902,pc202,pchn2,pchnf2,pcbn2,pcbnf2)

  }
  save(DATA,file = paste0("res-all/data-set-",i,"-S",S,".RData"))
  
}}

parallel::stopCluster(cl)

