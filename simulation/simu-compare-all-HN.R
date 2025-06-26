##
## Compare all methods based on HN model data generating process
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
s = c(20,50) ## #of population studies
set = expand.grid(
  t.theta = c(-0.7), ## true theta
  t.tau = sqrt(c(0.01, 0.3, 0.7)), ## true tau  
  t.rho = c(0.9), ## for population data rho does not matter the estimates
  n.median = c(100), ## median number of total subjects,
  grp.r = c(1), ##, 2group ratio: treat:control
  pmax = 0.99,
  pmin = 0.1
) %>% arrange(t.theta,n.median)
# set$p0 = 0.05
# set$t.rho = ifelse(set$t.theta>0, 0.8, -0.8)
# set$ymin = ifelse(set$n.median==50, 0, 0)
# set$ymax = ifelse(set$n.median==50, 5, 10)
set$ymin = 0
set$ymax = 10
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
  DATA = foreach(r=1:1000, .combine = rbind,.packages=c("tidyr","mnormt","dplyr","metafor"))  %dorng%  {

    set.gr = set[i,]
    pmax = set.gr$pmax
    pmin = set.gr$pmin
    
    ## create datasets
    plist = gen.data1(
      s=S, 
      n.med=set.gr$n.median,
      # n_min = 10, n_max=20,
      gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      y_min=set.gr$ymin,y_max=set.gr$ymax,
      Pnmax = pmax, Pnmin = pmin)

    
    ## selected data and population data
    sdata = plist$s.dt
    pdata = plist$p.dt

    ## data with n and lnOR
    data = lapply(list(pdata, sdata),
                  function(data) escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data)
    )
    lpdata = data[[1]]
    lsdata = data[[2]]
    
    # n.00 = lsdata %>%filter(y1 == 0, y0 == 0) %>%nrow()
    p.small = lpdata %>%summarise(prop = mean(y1 <=3 & y0 <=3))%>%c()
    s.small = lsdata %>%summarise(prop = mean(y1 <=3 & y0 <=3))%>%c()

    ## population and selected models -------
    
    ## set parset list
    parset.nn = list(
      mu.bound = abs(set.gr$t.theta)*1.5, 
      tau.bound = 1,
      eps = 1e-4,
      init.vals = c(set.gr$t.theta/2,set.gr$t.tau/2)
    )
    
    ## estimation without/with PB: NN, HN-GLMM, BN-GLMM on pdata and sdata
    fit.nn = lapply(
      list(lpdata,lsdata), 
      function(data) with(data, NN_LMM(yi, vi, parset=parset.nn)))
    
    pnn = c(fit.nn[[1]]$mu, fit.nn[[1]]$tau, rho=NA, rho.se=NA,
             cv = ifelse(is.null(fit.nn[[1]]$opt$convergence), NA, fit.nn[[1]]$opt$convergence))
    snn = c(fit.nn[[2]]$mu, fit.nn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[2]]$opt$convergence), NA, fit.nn[[2]]$opt$convergence))
    
    ## initial values for GLMM
    parset.glmm = list(
      mu.bound = abs(set.gr$t.theta)*1.5,
      tau.bound = 1,
      eps = 1e-4,
      integ.limit = 5, 
      cub.tol = 1e-3,
      init.vals = c(set.gr$t.theta/2,set.gr$t.tau/2)
    )
    
    fit.hn = lapply(
      list(lpdata, lsdata), 
      function(data) with(data, HN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    phn = c(fit.hn[[1]]$mu, fit.hn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[1]]$opt$convergence))
    shn = c(fit.hn[[2]]$mu, fit.hn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[2]]$opt$convergence), NA, fit.hn[[2]]$opt$convergence))
    
    fit.bn = lapply(
      list(lpdata, lsdata), 
      function(data) with(data, BN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    pbn = c(fit.bn[[1]]$mu, fit.bn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[1]]$opt$convergence), NA, fit.bn[[1]]$opt$convergence))
    sbn = c(fit.bn[[2]]$mu, fit.bn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[2]]$opt$convergence), NA, fit.bn[[2]]$opt$convergence))
   
    
    ## adjusted models -------
    ## Copas methods initial values
    parset.copas = list(
    mu.bound = abs(set.gr$t.theta)*1.5,
    tau.bound = 1,
    estimate.rho = TRUE, 
    eps = 1e-4,
    init.vals = c(set.gr$t.theta/2,set.gr$t.tau/2,set.gr$t.rho/2) ## initials for mu tau and rho
    )
    ## proposed method initial values
    parset.new = list(
      mu.bound = abs(set.gr$t.theta)*1.5,
      tau.bound = 1,
      estimate.rho = TRUE, 
      eps = 1e-4,
      integ.limit = 5, 
      cub.tol = 1e-3,
      init.vals = c(set.gr$t.theta/2,set.gr$t.tau/2,set.gr$t.rho/2)
    )
    
    nmin = min(lsdata$n)
    nmax = max(lsdata$n)
    
    
    func_list = list(
      function(x) with(x, 
                       COPAS1999(yi, n, Pnmax = pmax, Pnmin = pmin,
                                 parset=parset.copas)),
      function(x) with(x, 
                       COPAS2000(yi, vi, Psemax = pmax, Psemin = pmin,
                                 parset=parset.copas)),
      function(x) with(x,
                       COPAS_HNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                    n_min = nmin, n_max = nmax,
                                    parset=parset.new)),
      function(x) with(x,
                       COPAS_BNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                    n_min = nmin, n_max = nmax,
                                    parset=parset.new))
      
    )
    
    adj.list = lapply(func_list, function(f) f(lsdata))
    
    copas99 = c(adj.list[[1]]$mu, adj.list[[1]]$tau,adj.list[[1]]$rho,
              cv = ifelse(is.null(adj.list[[1]]$opt$convergence), NA, adj.list[[1]]$opt$convergence))
    copas20 = c(adj.list[[2]]$mu, adj.list[[2]]$tau,adj.list[[2]]$rho,
                cv = ifelse(is.null(adj.list[[2]]$opt$convergence), NA, adj.list[[2]]$opt$convergence))
    adjhn   = c(adj.list[[3]]$mu, adj.list[[3]]$tau,adj.list[[3]]$rho,
                cv = ifelse(is.null(adj.list[[3]]$opt$convergence), NA, adj.list[[3]]$opt$convergence))
    adjbn   = c(adj.list[[4]]$mu, adj.list[[4]]$tau,adj.list[[4]]$rho,
                cv = ifelse(is.null(adj.list[[4]]$opt$convergence), NA, adj.list[[4]]$opt$convergence))
  
    
    res = cbind(rbind(pnn,phn,pbn,snn,shn,sbn,copas99,copas20,adjhn,adjbn),
                p.prop=c(p.small, rep(NA,9)),
                s.prop=c(s.small, rep(NA,9)))
    res

  }
  save(DATA,file = paste0("res-all-HN/data-set-",i,"-S",S,".RData"))
  
}}

parallel::stopCluster(cl)

