#
library(dplyr)
library(MASS)
library(mnormt)
library(ggplot2)
library(foreach)
library(doRNG)
library(metafor)

s = c(10, 50) ## #of population studies

par(mfrow=c(2,2))
## summary1
sum.res.all = NULL
for(S in s[1]){
  for(i in 1:6){
    load(paste0("res-pop/data-set-",i,"-S",S,".RData"))

    DATA1=DATA%>%t()
    dim(DATA1) = c(6,6,1000)

    ## remove nonconverged values
    for(i in 1:1000){
      DATA1[,,i][1:5,(DATA1[,,i][6,]!=0)]=NA
    }

    sum.res = apply(DATA1, c(1, 2), function(x) median(x, na.rm=T))
    colnames(sum.res) = c("pnn1","pnn2","phn1","phn2","pbn1","pbn2")
    rownames(sum.res) = c("mu","mu.se","tau","tau.se","tau2","cv")
    sum.res.all=rbind(sum.res.all, sum.res)

  }}

## data generating 1
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 6), c(1,3,5)]-(-0.7)
matplot(df.mu, type="b")
abline(h=0)
title("S=10 and data generating 1")

## data generating 2
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 6), c(2,4,6)]-(-0.7)
matplot(df.mu, type="b")
abline(h=0)
title("S=10 and data generating 2")


## summary2
sum.res.all = NULL
for(S in s[2]){
  for(i in 1:6){
    load(paste0("res-pop/data-set-",i,"-S",S,".RData"))

    DATA1=DATA%>%t()
    dim(DATA1) = c(6,6,1000)

    ## remove nonconverged values
    for(i in 1:1000){
      DATA1[,,i][1:5,(DATA1[,,i][6,]!=0)]=NA
    }

    sum.res = apply(DATA1, c(1, 2), function(x) median(x, na.rm=T))
    colnames(sum.res) = c("pnn1","pnn2","phn1","phn2","pbn1","pbn2")
    rownames(sum.res) = c("mu","mu.se","tau","tau.se","tau2","cv")
    sum.res.all=rbind(sum.res.all, sum.res)

  }}

## data generating 1
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 6), c(1,3,5)]-(-0.7)
matplot(df.mu, type="b")
abline(h=0)
title("S=50 and data generating 1")

## data generating 2
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 6), c(2,4,6)]-(-0.7)
matplot(df.mu, type="b")
abline(h=0)
title("S=50 and data generating 2")

par(mfrow=c(1,1))