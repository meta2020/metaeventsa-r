library(data.table)

data.gen = function(param){
  M = param [1]
  tlogor = param [3]
  ttau = param [2]
  alpha = param [4]
  beta = param [5]
  
  rbindlist(lapply(1:1000,function(i){
    tlogor.i  = rnorm( M, tlogor, ttau ) # empirical logORi
    tOR.i     = exp( tlogor.i ) # empirical odds ratio
    #----------generate a b c d of OR,then calculate se----------------#
    Pc                   = runif(M,0.2,0.9) # Probability of event in the control group
    m                    = tOR.i * Pc / ( 1 - Pc )
    Pt                   = m / ( 1 + m )   # calculation of probability of event in the experimental group
    u                    = 2*round(rlnorm(M,5,1)/2,0) # total sample size
    nt  = nc = ifelse(u<20,10,u/2)
    all.dat = rbindlist(lapply(1:M,function(i){
      a = sum(rbinom(nt[i],1,Pt[i]))
      b = nt[i] -a
      c = sum(rbinom(nc[i],1,Pc[i]))
      d = nc[i]-c
      if ( 0 %in% c( a,b,c,d )) {
        a = a + 0.5
        b = b + 0.5
        c = c + 0.5
        d = d + 0.5
      }
      logOR = log(a*d/(b*c))
      se = sqrt(1/a+1/b+1/c+1/d)
      ni = nt[i]+nc[i]
      data.table(logOR,se,ni)}))
    all.dat[,wi:=.(pnorm(alpha+beta*logOR/se))]
    all.dat[,zi:=.(rbinom (M,1,wi))]
    data.table(all.dat[,.(logOR=logOR,se=se,zi=zi,ni=ni)],t(param),i)}))
}

