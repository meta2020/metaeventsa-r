#'
#' Meta-analysis of the logit-proportion (Example 2)
#'
#' Load R functions
rm(list=ls())
file.sources = list.files("Rfunc/")
sapply(paste0("Rfunc/", file.sources), source)

#' Data are from 1. Stijnen T, Hamza TH, Özdemir P. 
#' Random effects meta-analysis of event outcome in the framework of the 
#' generalized linear mixed model with applications in sparse data. 
#' Stat Med. 2010;29(29):3046-3067. 
data = read.csv("niel-weise21.csv")

#' Derive continuous outcomes (logit-proportion and se)
metadt = meta::metaprop(event=y1, n=n1, data=data, 
                  sm="PLOGIT", method.incr = "only0", allstudies = T)
data$yi = metadt$TE
data$vi = (metadt$seTE)^2

#' Meta-analysis without PB ----------
#' Data
y1 = data$y1
n1 = data$n1
yi = data$yi
vi = data$vi

#' NN model
lgtP_nn = NN_LMM(
  yi=yi, vi=vi, 
  parset=list(
    mu.bound = 10, 
    tau.bound = 5,
    eps = 1e-3,
    init.vals = NULL
  ))
lgtP_nn_lb = lgtP_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lgtP_nn$mu[2]
lgtP_nn_ub = lgtP_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lgtP_nn$mu[2]
sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lgtP_nn$mu[1], lgtP_nn_lb, lgtP_nn_ub, lgtP_nn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lgtP_nn$tau[1], lgtP_nn$tau[2])

#' BN-GLMM
lgtP_bn = BN_GLMM_prop(
  y1 = y1, n1 = n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lgtP_bn_lb = lgtP_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lgtP_bn$mu[2]
lgtP_bn_ub = lgtP_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lgtP_bn$mu[2]
sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lgtP_bn$mu[1], lgtP_bn_lb, lgtP_bn_ub, lgtP_bn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lgtP_bn$tau[1], lgtP_bn$tau[2])


#' Proposed sensitivity analysis for PB in meta-analysis ----------
#' (Pnmax = 0.999, Pnmin = p_sa)
#' (Psemax = p_sa, Psemin = 0.999)
#' 
#' Proposal for BN-GLMM
p_sa = c(0.99, seq(0.9, 0.1, -0.1))
lgtP_COPAS_BNGLMM = vapply(
  p_sa, 
  function(p) {
    mod = COPAS_BNGLMM_prop(
      y1=y1, n1=n1, Pnmax = 0.999, Pnmin = p, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        integ.limit = 10, 
        cub.tol = 1e-5,
        init.vals = NULL)
    )
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$a, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "a0"=0, "a1"=0, "converge"=0))
colnames(lgtP_COPAS_BNGLMM) = paste0("p = ", p_sa)
lgtP_COPAS_BNGLMM

#' Copas-Heckman-type selection model (2000)
lgtP_COPAS2000 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS2000(
      yi = yi, vi = vi, Psemax = p, Psemin = 0.999, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS2000) = paste0("p = ", p_sa)
lgtP_COPAS2000

#' Copas-N-type selection model (2000)
lgtP_COPAS1999 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS1999(
      yi = yi, ni = n1, Pnmax = 0.999, Pnmin = p, 
      parset = list(
          mu.bound = 10,
          tau.bound = 5,
          estimate.rho = TRUE, 
          eps = 1e-3,
          init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS1999) = paste0("p = ", p_sa)
lgtP_COPAS1999

#' NUMBER OF THE UNPUBLISHED
M_propos = sapply(p_sa, function(p) {
  
  P_max = 0.999
  P_min = p

  n_min = min(n1) 
  n_max = max(n1)
  
  a1 = (qnorm(P_max)-qnorm(P_min))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(P_max)-a1*sqrt(n_max)
  sum((1 - pnorm(a0+a1*sqrt(n1)))/pnorm(a0+a1*sqrt(n1))) 
  
})


M_copas = sapply(p_sa, function(p) {
  
  P_max = p
  P_min = 0.999

  se_min_inv = 1/sqrt(min(vi)) 
  se_max_inv = 1/sqrt(max(vi))
  
  gamma1 = (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 = qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi)))/pnorm(gamma0+gamma1/sqrt(vi)))
  
})


#' Table 1
tab1 = data.frame(
  BN = t(lgtP_COPAS_BNGLMM[c(1,3,4),]),
  CH = t(lgtP_COPAS2000[c(1,3,4),]), # copas-heckman
  CN = t(lgtP_COPAS1999[c(1,3,4),]), # copas-n
  M.c = round(M_copas), 
  M.p = round(M_propos)
)

tab1_p1 = c(
  BN = lgtP_bn$mu[1], BN.mu.lb = unname(lgtP_bn_lb), BN.mu.ub = unname(lgtP_bn_ub),
  CH = lgtP_nn$mu[1], NN.mu.lb = unname(lgtP_nn_lb), NN.mu.ub = unname(lgtP_nn_ub),
  CN = lgtP_nn$mu[1], NN.mu.lb = unname(lgtP_nn_lb), NN.mu.ub = unname(lgtP_nn_ub),
  M.c = 0, M.p = 0
)

# tab1_all = rbind("p = 1" = tab1_p1, tab1)
tab1_all = tab1
tab1_all$pnmin = p_sa
tab1_all$pnmax = rep(0.999, 10)

#' SAVE RESULTS1
#' 
save(lgtP_COPAS_BNGLMM, lgtP_COPAS2000,lgtP_COPAS1999,
     M_propos, M_copas, tab1_all,
     lgtP_bn, lgtP_nn,
     file = "example-supp-bias.RData")


