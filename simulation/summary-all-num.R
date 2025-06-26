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

## summary1
sum.res.all = NULL
S=10
for(i in 1:6){
  load(paste0("res-all/data-set-",i,"-S",S,".RData"))
  DATA0 = DATA%>%as.data.frame()
  DATA0$gen = rep(rep(1:2,each=12),rtimes)
  
  DATA1=DATA0[DATA0$gen==1,1:8]%>%t()
  dim(DATA1) = c(8,12,rtimes)
  
  ## remove nonconverged values
  for(j in 1:rtimes){
    DATA1[1:7,,j][,(DATA1[,,j][8,]!=0)]=NA
  }
  
  sum.res = apply(DATA1, c(1, 2), function(x) median(x, na.rm=T))
  colnames(sum.res) = c("pnn1","phn1","pbn1","snn1","shn1","sbn1","pc901","pc201","pchn1","pchnf1","pcbn1","pcbnf1")
  rownames(sum.res) = c("mu","mu.se","tau","tau.se","tau2","rho","rho.se","cv")
  sum.res.all=rbind(sum.res.all, sum.res)
}

#bias on mu
df.mu = sum.res.all[seq(1, nrow(sum.res.all), by = 8), ]-t.mu
# Select columns with names containing "1"
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("nn","0")))
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("hn")))
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("bn")))
