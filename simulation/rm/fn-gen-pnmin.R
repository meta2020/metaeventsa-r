#' Estimate Pmin in each meta-analysis
#'
#' @param s number of studies
#' @param theta true log OR
#' @param tau true between-study standard deviation
#' @param rho true average control group probability
#' @param n_min minimum value for total subjects
#' @param n_max maximum value for total subjects
#' @param y_min minimum value for events (ycontrol+ytreat)
#' @param y_max maximum value for events (ycontrol+ytreat)
#' @param Pnmax probabilities of publishing a study with maximum subjects
#' @param p overall selection probbaility
#' 
gen.pnmin1 = function(s,theta,tau,rho,
                    n_min=200,n_max=400,
                    y_min=10,y_max=20,
                    Pnmax = 0.99,p=0.7){
  
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  n0i = runif( s, min = n_min, max = n_max ) %>% round()
  n1i = runif( s, min = n_min, max = n_max ) %>% round()
  yi  = runif( s, min = y_min, max = y_max ) %>% round()
  study.pi = plogis( log(n1i/n0i) + thetai )
  
  y1i = rbinom( s, size = yi, prob = study.pi ) 
  y0i = yi- y1i 
  
  p.dt = data.frame(y1=y1i,y0=y0i,n1=n1i,n0=n0i)
  
  ## SELECTION
  ni=n0i+n1i
  
  n_min=min(ni) 
  n_max=max(ni)
  
  est.pnmin=function(Pnmin){
    
  a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0=qnorm(Pnmax)-a1*sqrt(n_max)
  
  zi=a0+a1*sqrt(ni)+deltai
  p.dt$z=1*(zi>0)
  
  pz=pnorm(a0+a1*sqrt(ni))
  
  return(mean(pz,na.rm=T)-p)
  
  }
  
  pnmin=uniroot(est.pnmin, interval = c(0.001,0.999), extendInt="no")$root
  
  return(pnmin)
  
}

# gen.pnmin(s=5000,theta=0.4,tau=0.15,rho=0.8,
#           n_min=10,n_max=20,y_min=3,y_max=5,
#           Pnmax = 0.99,p=0.7)


#' @param s number of studies
#' @param theta true log OR
#' @param tau true between-study standard deviation
#' @param rho true average control group probability
#' @param p0i.range range of probability of event in the control group
#' @param n.median median of sample size in the control group
#' @param gr group ratio (ntreat/ncontrol): 1,2
#' @param Pnmax probabilities of publishing a study with maximum subjects
#' @param Pnmin probabilities of publishing a study with minimum subjects

gen.pnmin2 <- function(s,theta,tau,rho,
                      p0i.range,
                      n.median,gr,
                      Pnmax, p) {
  
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  thetai = rnorm( s, theta, tau ) # empirical logORi
  ORi    = exp( thetai ) # empirical odds ratio
  
  p0i = runif(s, p0i.range[1], p0i.range[2]) # Probability of event in the control group
  p1i  = p0i*ORi / ( 1 -p0i + p0i*ORi )   # calculation of probability of event in the experimental group
  
  n0i = round(rlnorm(s,log(n.median),1)) 
  n0i = ifelse(n0i<10,10,n0i)
  n1i = n0i*gr
  
  p.dt = vapply(1:s,function(i){
    y1i = sum(rbinom(n1i[i],1,p1i[i]))
    y0i = sum(rbinom(n0i[i],1,p0i[i]))
    n1i = n1i[i]
    n0i = n0i[i]
    c(y1i, y0i, n1i, n0i)
  },c("y1"=0,"y0"=0,"n1"=0,"n0"=0))%>%t()%>%as.data.frame()  ## datasets
  
  ## selective publication
  ni=n0i+n1i
  
  n_min=min(ni) 
  n_max=max(ni)
  
  est.pnmin=function(Pnmin){
    
    a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
    a0=qnorm(Pnmax)-a1*sqrt(n_max)
    
    zi=a0+a1*sqrt(ni)+deltai
    p.dt$z=1*(zi>0)
    
    p.dt$pz=pnorm(a0+a1*sqrt(ni))
    
    return(mean(p.dt$pz,na.rm=T)-p)
    
  }
  
  pnmin=uniroot(est.pnmin, interval = c(0.001,0.999), extendInt="no")$root
  
  return(pnmin)
  

  
  p.dt$pz= pnorm(a0+a1*sqrt(ni))
  p=mean(p.dt$pz)
  
  s.dt = p.dt[p.dt$z>0,]
  M  = sum((1-s.dt$pz)/s.dt$pz)%>%round()  ## theoretical # of missing studies 
  M.e = s-sum(p.dt$z) ## empirical # of missing studies 
  
  return(list(p.dt=p.dt,s.dt=s.dt,
              m=M, e.m=M.e,
              a1=a1,a0=a0,p=p))
}

gdt=gen.data2(s=20,theta=1.5,tau=0.5,rho=0.8,
              p0i.range=c(0.1,0.2),
              n.median=10,gr=1,
              Pnmax=0.99, Pnmin=0.5)
