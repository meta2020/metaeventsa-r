# s=10
# n.med=50
# gr=1
# theta=-2
# tau=0.05
# rho=0.8
# y_min=0
# y_max=10
# Pnmax = 0.9999
# Pnmin = 0.2

gen.data2 = function(
    s, 
    # n_min, n_max,
    n.med,
    gr,
    theta,tau,rho,
    y_min,y_max,
    Pnmax, Pnmin){
  
  # n0i = runif( s, min = n_min, max = n_max ) %>% round()
  # n1i = n0i*gr
  # n = runif( s, min = n_min, max = n_max )
  n = round(rlnorm(s,log(n.med),1)) # total #subjects
  n = ifelse(n<20,20,n)
  
  # n0i = rbinom(s, n, 1 / (1 + gr))%>% round() # #control subjects
  
  n0i = (n*(1 / (1 + gr))) %>% round()
  n1i = n-n0i # #treatment subjects
  ni=n0i+n1i
  
  # generate deltai and thetai
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  # generate yi
  yi  = runif(s, min =y_min, max = y_max) %>% round()
  yi[which(yi>ni)]=ni[which(yi>ni)]
  # study.pi = plogis(log(n1i/n0i) + thetai)
  
  # y1i = rbinom(s, size = yi, prob = study.pi)
  # yi = BiasedUrn::rFNCHypergeo(nran=s, m1=n1[i], m2=n0[i], n=yi[i], odds=exp(thetai))
  # y1i = sapply(1:s, function(i) MCMCpack::rnoncenhypergeom(1, n1i[i], n0i[i], yi[i], exp(thetai[i])))
  y1i = sapply(1:s, function(i) rbinom(1, yi[i], plogis(log(n1i[i]/n0i[i])+thetai[i])) )
  
  y0i = yi- y1i 
  y1i[which(y1i>n1i)]=n1i[which(y1i>n1i)]
  y0i[which(y0i>n0i)]=n0i[which(y0i>n0i)]
  
  p.dt = data.frame(y1=y1i,y0=y0i,n1=n1i,n0=n0i,n=ni)
  
  ## selective process
  n_min=min(ni) 
  n_max=max(ni)
  
  a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0=qnorm(Pnmax)-a1*sqrt(n_max)
  
  zi=a0+a1*sqrt(ni)+deltai
  z=1*(zi>0)
  
  pz=pnorm(a0+a1*sqrt(ni))
  # z=sapply(1:s, function(i) rbinom(1,1,pz[i]))
  pz.s = pz[z>0]
  M= sum((1-pz.s)/pz.s)%>%round()
  # p.dt$z1=z1
  # p.dt$zi=zi
  p.dt$z=z
  p.dt$pz=pz
  
  s.dt=p.dt[z>0,]
  
  res.list = list(
    p.dt=p.dt,
    s.dt=s.dt,
    a=c(a1=a1,a0=a0),
    Miss=c(M.e=M, M.o=sum(z==0)), 
    p=c(p.e=mean(pz), p.o=mean(z>0))
  )
  return(res.list)
  
}