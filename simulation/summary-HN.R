library(dplyr)
library(MASS)
library(mnormt)
library(ggplot2)
library(foreach)
library(doRNG)
library(metafor)

# s = c(10, 50) ## #of population studies
rtimes=200
t.mu = -2

## summary1
sum.res.all = NULL
sum.cv.all = NULL
S=50
for(i in 1:3){
  
  load(paste0("res-all-HN//data-set-",i,"-S",S,".RData"))
  DATA0 = DATA %>% t()%>% as.numeric() %>% 
    array(., dim = c(10, 10, 200),
          dimnames = list(colnames(DATA),rownames(DATA)[1:10],c(1:200)))
  # remove nonconverged values
  for(j in 1:rtimes){
    DATA0[1:7,,j][,(DATA0[,,j][8,]!=0)]=NA
  }
  mean_matrix <- apply(DATA0, c(1, 2), function(x) mean(x, na.rm=T))[8,]
  median_matrix <- apply(DATA0, c(1, 2), function(x) median(x, na.rm=T))[1,]
  sum.res.all=rbind(sum.res.all, median_matrix)
  sum.cv.all=rbind(sum.cv.all, mean_matrix)
}

#bias on mu
df.mu = sum.res.all-t.mu

# Select columns with names containing "1"
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("nn","0")))
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("hn")))
df.mu%>%as.data.frame() %>% dplyr::select(matches(c("bn")))
