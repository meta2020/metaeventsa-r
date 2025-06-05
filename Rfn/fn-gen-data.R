#' Generate meta-analysis with sparse event
#'
#' generate data from BN model
#' @param s number of studies
#' @param theta true log OR
#' @param tau true between-study standard deviation
#' @param rho true average control group probability
#' @param n_min minimum value for total subjects
#' @param n_max maximum value for total subjects
#' @param y_min minimum value for events (ycontrol+ytreat)
#' @param y_max maximum value for events (ycontrol+ytreat)
#' @param Pnmax probabilities of publishing a study with maximum subjects
#' @param Pnmin probabilities of publishing a study with minimum subjects
#' 
#' s=10
#' n.med=10
#' gr=1
#' tau=0.5
#' rho=0.8
#' theta=-0.7
gen.data1 = function(
  s, n.med, gr,
  theta,tau,rho,
  y_min=10,y_max=20,
  Pnmax = 0.99, Pnmin = 0.5){

  # n0i = runif( s, min = n_min, max = n_max ) %>% round()
  # n1i = runif( s, min = n_min, max = n_max ) %>% round()
  
  n = round(rlnorm(s,log(n.med),1)) # total #subjects
  n = ifelse(n<10,10,n)
  n0i = rbinom(s, n, 1 / (1 + gr)) # #control subjects
  n1i = n - n0i # #treatment subjects

  # generate deltai and thetai
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]

  # generate yi
  yi  = runif(s, min = y_min, max = y_max) %>% round()
  study.pi = plogis(log(n1i/n0i) + thetai)
  
  y1i = rbinom(s, size = yi, prob = study.pi) 
  y0i = yi- y1i 
  ni=n0i+n1i
  p.dt = data.frame(y1=y1i,y0=y0i,n1=n1i,n0=n0i,n=ni)
  
  ## selective process
  n_min=min(ni) 
  n_max=max(ni)
  
  a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0=qnorm(Pnmax)-a1*sqrt(n_max)
  
  zi=a0+a1*sqrt(ni)+deltai
  p.dt$z=1*(zi>0)
  
  pz=pnorm(a0+a1*sqrt(ni))
  p.dt$pz=pz
  p=mean(pz)
  
  s.dt = p.dt[p.dt$z>0,]
  M=sum((1-pz)/pz)%>%round()
  M.e = s-sum(p.dt$z)
  
  return(list(
    p.dt=p.dt,s.dt=s.dt,
    m=M, e.m=M.e,
    a1=a1,a0=a0,p=p))
  
}


# 
#' Generate data refering to Jansen K, Holling H. 
#' Random-effects meta-analysis models for the odds ratio in the case of 
#' rare events under different data-generating models: A simulation study. 
#' Biometrical J. 2023;65(3):1-28. 
#'
#' @param s number of studies
#' @param theta true log OR
#' @param tau true between-study standard deviation
#' @param rho true average control group probability
#' @param p0i.range range of probability of event in the control group
#' @param n0.median median of sample size in the control group
#' @param gr group ratio (ntreat/ncontrol): 1,2
#' @param Pnmax probabilities of publishing a study with maximum subjects
#' @param Pnmin probabilities of publishing a study with minimum subjects
#' 
gen.data2 = function(
    s,theta,tau,rho,
    p0,n.med,gr,
    Pnmax, Pnmin) {

  ## number of subjects
  n = round(rlnorm(s,log(n.med),1))
  n = ifelse(n<10,10,n)
  n0i = rbinom(s, n, 1 / (1 + gr)) # decides for each observation
  n1i = n - n0i
  
  ## thetai=log(ORi) and deltai 
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  ## empirical ORi
  ORi = exp(thetai)

  ## true p0
  mu0 = qlogis(p0) ## logit(p0) baseline
  
  ## empirical logit(p0i)
  logitp0i = rnorm(s, mu0, tau^2 / sqrt(2))
  # logitp1i = rnorm(s, mu0+theta, tau^2 / sqrt(2))
  
  ## empirical logit(p1i): thetai = logitp1i-logitp0i
  logitp1i = thetai+logitp0i

  p0i = plogis(logitp0i) ## logistic(logitp0i)
  p1i = plogis(logitp1i)
  
  ## number of events
  y0i = rbinom(s, n0i, p0i)
  y1i = rbinom(s, n1i, p1i)

  p.dt = data.frame(y0=y0i,y1=y1i,n0=n0i,n1=n1i,n=n)
  
  ## selective publication
  n_min=min(n) 
  n_max=max(n)
  a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0=qnorm(Pnmax)-a1*sqrt(n_max)
  
  zi=a0+a1*sqrt(n)+deltai
  p.dt$z=1*(zi>0)
  
  p.dt$pz= pnorm(a0+a1*sqrt(n))
  p=mean(p.dt$pz)
  
  s.dt = p.dt[p.dt$z>0,]
  M  = sum((1-s.dt$pz)/s.dt$pz)%>%round()  ## theoretical # of missing studies 
  M.e = s-sum(p.dt$z) ## empirical # of missing studies 
  # 
  # n_min=min(s.dt$n) 
  # n_max=max(s.dt$n)
  # sa1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  # sa0=qnorm(Pnmax)-a1*sqrt(n_max)

  return(list(p.dt=p.dt,s.dt=s.dt,
              m=M, e.m=M.e,
              a1=a1,a0=a0,p=p))
}


##
## From Ao
##

gen.data3(
    s=10, n.med=10, gr=1,
    theta=-0.7,tau=0.01,rho=0.8,
    p0=0.05,
    Pnmax=0.99, Pnmin=0.5
    )
gen.data3 = function(
    s, n.med, gr,
    theta,tau,rho,
    p0,
    Pnmax, Pnmin){
  
    n = round(rlnorm(s,log(n.med),1)) # total #subjects
    n = ifelse(n<10,10,n)
    n0i = rbinom(s, n, 1 / (1 + gr)) # #control subjects
    n1i = n - n0i # #treatment subjects
    
    ## thetai=log(ORi) and deltai 
    sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
    m = MASS::mvrnorm(s,c(theta,0),sigma)
    thetai=m[,1]
    deltai=m[,2]
  
    tOR.i = exp(thetai) # empirical odds ratio
    #----------generate a b c d of OR,then calculate se----------------#
    # Pc                   = runif(M,0.2,0.9) # Probability of event in the control group
    Pci = plogis(rnorm(s, qlogis(p0), tau^2 / sqrt(2)))
    m                    = tOR.i * Pci / ( 1 - Pci )
    Pti                   = m / ( 1 + m )   # calculation of probability of event in the experimental group
    # u                    = 2*round(rlnorm(s,5,1)/2,0) # total sample size
    u = round(rlnorm(s,log(n.med),1))
    nc= rbinom(s, u, 1 / (1 + gr))
    nt= u-nc
    # nt  = ifelse(u<20,10, 1/(1 + gr))

    all.dat = vapply(1:s,function(i){
      a = sum(rbinom(nt[i],1,Pti[i]))
      b = nt[i]
      c = sum(rbinom(nc[i],1,Pci[i]))
      d = nc[i]
      c(a,b,c,d)
      },c(y1i=0,n1i=0,y0i=0,n0i=0))%>%t()%>%data.frame()
    
    ni = all.dat$n1i+all.dat$n0i
    p.dt = data.frame(y1=all.dat$y1i,y0=all.dat$y0i,
                      n1=all.dat$n1i,n0=all.dat$n0i,n=ni)
    
    ## selective process
    n_min=min(ni) 
    n_max=max(ni)
    
    a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
    a0=qnorm(Pnmax)-a1*sqrt(n_max)
    
    zi=a0+a1*sqrt(ni)+deltai
    p.dt$z=1*(zi>0)
    
    pz=pnorm(a0+a1*sqrt(ni))
    p.dt$pz=pz
    p=mean(pz)
    
    s.dt = p.dt[p.dt$z>0,]
    M=sum((1-pz)/pz)%>%round()
    M.e = s-sum(p.dt$z)
    
    return(list(
      p.dt=p.dt,s.dt=s.dt,
      m=M, e.m=M.e,
      a1=a1,a0=a0,p=p))
    
    
}

