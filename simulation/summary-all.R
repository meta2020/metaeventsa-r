#
library(dplyr)
library(MASS)
library(mnormt)
library(ggplot2)
library(foreach)
library(doRNG)
library(metafor)

s = c(10, 50) ## #of population studies
rtimes=200
t.mu = -0.7
par(mfrow=c(2,2))
## summary1
sum.res.all = NULL
for(S in s[1]){
  for(i in 1:6){
    load(paste0("res-all/data-set-",i,"-S",S,".RData"))
    DATA0 = DATA%>%as.data.frame()
    DATA0$gen = rep(rep(1:2,each=12),rtimes)
    
    DATA1=DATA0[DATA0$gen==2,1:8]%>%t()
    dim(DATA1) = c(8,12,rtimes)

    ## remove nonconverged values
    for(j in 1:rtimes){
      DATA1[,,j][,(DATA1[,,j][8,]!=0)]=NA
    }

    sum.res = apply(DATA1, c(1, 2), function(x) median(x, na.rm=T))
    colnames(sum.res) = c("pnn1","phn1","pbn1","snn1","shn1","sbn1","pc901","pc201","pchn1","pchnf1","pcbn1","pcbnf1")
    rownames(sum.res) = c("mu","mu.se","tau","tau.se","tau2","rho","rho.se","cv")
    sum.res.all=rbind(sum.res.all, sum.res)

  }}

## data generating 1
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 8), ]-t.mu
matplot(df.mu[,7:12], type="b")
abline(h=0)
title("S=10 and data generating 1")

## data generating 2
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 8), c(2,4,6,8,10,12)]-t.mu
matplot(df.mu[,c(1,2,3,5)], type="b")
abline(h=0)
title("S=10 and data generating 2")


## summary2
sum.res.all = NULL
for(S in s[2]){
  for(i in 1:6){
    load(paste0("res-adjS/data-set-",i,"-S",S,".RData"))
    
    DATA1=DATA%>%t()
    dim(DATA1) = c(8,12,rtimes)
    
    ## remove nonconverged values
    for(j in 1:rtimes){
      DATA1[,,j][,(DATA1[,,j][8,]!=0)]=NA
    }
    
    sum.res = apply(DATA1, c(1, 2), function(x) median(x, na.rm=T))
    colnames(sum.res) = c("pc901","pc902",
                          "pc201","pc202",
                          "pchn1","pchn2",
                          "pchnf1","pchnf2",
                          "pcbn1","pcbn2",
                          "pcbnf1","pcbnf2")
    rownames(sum.res) = c("mu","mu.se","tau","tau.se","tau2","rho","rho.se","cv")
    sum.res.all=rbind(sum.res.all, sum.res)
    
  }}

## data generating 1
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 8), c(1,3,5,7,9,11)]-t.mu
matplot(df.mu[,c(1,2,3,5)], type="b")
abline(h=0)
title("S=50 and data generating 1")

## data generating 2
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 8), c(2,4,6,8,10,12)]-t.mu
matplot(df.mu[,c(1,2,3,5)], type="b")
abline(h=0)
title("S=50 and data generating 2")
par(mfrow=c(1,1))
