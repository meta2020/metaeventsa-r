
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
# gen.data2 = function(
#     s,theta,tau,rho,
#     p0,n.med,gr,
#     Pnmax, Pnmin) {
# 
#   ## number of subjects
#   n = round(rlnorm(s,log(n.med),1))
#   n = ifelse(n<10,10,n)
#   n0i = rbinom(s, n, 1 / (1 + gr)) # decides for each observation
#   n1i = n - n0i
#   
#   ## thetai=log(ORi) and deltai 
#   sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
#   m = MASS::mvrnorm(s,c(theta,0),sigma)
#   thetai=m[,1]
#   deltai=m[,2]
#   
#   ## empirical ORi
#   ORi = exp(thetai)
# 
#   ## true p0
#   mu0 = qlogis(p0) ## logit(p0) baseline
#   
#   ## empirical logit(p0i)
#   logitp0i = rnorm(s, mu0, tau^2 / sqrt(2))
#   # logitp1i = rnorm(s, mu0+theta, tau^2 / sqrt(2))
#   
#   ## empirical logit(p1i): thetai = logitp1i-logitp0i
#   logitp1i = thetai+logitp0i
# 
#   p0i = plogis(logitp0i) ## logistic(logitp0i)
#   p1i = plogis(logitp1i)
#   
#   ## number of events
#   y0i = rbinom(s, n0i, p0i)
#   y1i = rbinom(s, n1i, p1i)
# 
#   p.dt = data.frame(y0=y0i,y1=y1i,n0=n0i,n1=n1i,n=n)
#   
#   ## selective publication
#   n_min=min(n) 
#   n_max=max(n)
#   a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
#   a0=qnorm(Pnmax)-a1*sqrt(n_max)
#   
#   zi=a0+a1*sqrt(n)+deltai
#   p.dt$z=1*(zi>0)
#   
#   p.dt$pz= pnorm(a0+a1*sqrt(n))
#   p=mean(p.dt$pz)
#   
#   s.dt = p.dt[p.dt$z>0,]
#   M  = sum((1-s.dt$pz)/s.dt$pz)%>%round()  ## theoretical # of missing studies 
#   M.e = s-sum(p.dt$z) ## empirical # of missing studies 
#   # 
#   # n_min=min(s.dt$n) 
#   # n_max=max(s.dt$n)
#   # sa1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
#   # sa0=qnorm(Pnmax)-a1*sqrt(n_max)
# 
#   return(list(p.dt=p.dt,s.dt=s.dt,
#               m=M, e.m=M.e,
#               a1=a1,a0=a0,p=p))
# }


##
## From Ao
##
gen.data3 = function(
    s, n.med, gr,
    theta,tau,rho,
    p0,
    Pnmax, Pnmin){
  
    n = round(rlnorm(s,log(n.med),1)) # total #subjects
    n = ifelse(n<10,10,n)
    n0 = rbinom(s, n, 1 / (1 + gr)) # #control subjects
    n1 = n - n0 # #treatment subjects
    
    ## thetai=log(ORi) and deltai 
    sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
    m = MASS::mvrnorm(s,c(theta,0),sigma)
    thetai=m[,1]
    deltai=m[,2]
  
    ORi = exp(thetai) # empirical odds ratio
    # pi0 = runif(s,0.2,0.9)
    pi0 = plogis(rnorm(s, qlogis(p0), tau / 2)) # prob of control in study i
    m = ORi * pi0 / ( 1 - pi0 )
    pi1 = m / ( 1 + m )   # prob of event in study i

    all.dat = vapply(1:s,function(i){
      yi1 = sum(rbinom(n1[i],1,pi1[i]))
      ni1 = n1[i]
      yi0 = sum(rbinom(n0[i],1,pi0[i]))
      ni0 = n0[i]
      c(yi1,ni1,yi0,ni0)
      },c(y1=0,n1=0,y0=0,n0=0))%>%t()%>%data.frame()
    
    n = all.dat$n1+all.dat$n0
    p.dt = data.frame(y1=all.dat$y1,y0=all.dat$y0,
                      n1=all.dat$n1,n0=all.dat$n0,n=n)
    
    ## selective process
    n_min=min(n) 
    n_max=max(n)
    
    a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
    a0=qnorm(Pnmax)-a1*sqrt(n_max)
    
    zi=a0+a1*sqrt(n)+deltai
    p.dt$z=1*(zi>0)
    
    pz=pnorm(a0+a1*sqrt(n))
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

