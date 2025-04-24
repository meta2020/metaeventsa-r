#'
#' Meta-analysis of the lnOR (Example 1)
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

#' Meta-analysis without PB ----------
#' Data
y1 = data$y1
y0 = data$y0
n1 = data$n1
n0 = data$n0
ni = n1+n0

#' Derive continuous outcomes (lnOR and se)
yvi = metafor::escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0)
yi = yvi[,1]
vi = yvi[,2]


#' NN model
lnOR_nn = NN_LMM(
  yi=yi, vi=vi, 
  parset=list(
    mu.bound = 10, 
    tau.bound = 5,
    eps = 1e-3,
    init.vals = NULL
  ))
lnOR_nn_lb = lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_nn$mu[2]
lnOR_nn_ub = lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_nn$mu[2]
sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_nn$mu[1], lnOR_nn_lb, lnOR_nn_ub, lnOR_nn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lnOR_nn$tau[1], lnOR_nn$tau[2])

#' HN-GLMM (takes longer time)
lnOR_hn = HN_GLMM(
  y0 = y0, y1 = y1, n0 = n0, n1 = n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lnOR_hn_lb = lnOR_hn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_hn$mu[2]
lnOR_hn_ub = lnOR_hn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_hn$mu[2]
sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_hn$mu[1], lnOR_hn_lb, lnOR_hn_ub, lnOR_hn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lnOR_hn$tau[1], lnOR_hn$tau[2])

#' BN-GLMM
lnOR_bn = BN_GLMM(
  y0 = y0, y1 = y1, n0 = n0, n1 = n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lnOR_bn_lb = lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_bn$mu[2]
lnOR_bn_ub = lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_bn$mu[2]
sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_bn$mu[1], lnOR_bn_lb, lnOR_bn_ub, lnOR_bn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lnOR_bn$tau[1], lnOR_bn$tau[2])



#' Proposed sensitivity analysis for PB in meta-analysis ----------
#' (Pnmax = 0.999, Pnmin = p_sa)
#' (Psemax = p_sa, Psemin = 0.999)
#' 
#' Proposal for HN-GLMM
p_sa = c(0.99, seq(0.9, 0.1, -0.1))
lnOR_COPAS_HNGLMM = vapply(
  p_sa, 
  function(p) {
    mod = COPAS_HNGLMM(
      y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.999, Pnmin = p, 
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
colnames(lnOR_COPAS_HNGLMM) = paste0("p = ", p_sa)
lnOR_COPAS_HNGLMM

#' Proposal for BN-GLMM
lnOR_COPAS_BNGLMM = vapply(
  p_sa, 
  function(p) {
    mod = COPAS_BNGLMM(
      y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.999, Pnmin = p, 
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
colnames(lnOR_COPAS_BNGLMM) = paste0("p = ", p_sa)
lnOR_COPAS_BNGLMM

#' Copas-Heckman-type selection model (2000)
lnOR_COPAS2000 = vapply(
  p_sa, 
  function(p) {
    mod = suppressWarnings(COPAS2000(
      yi = yi, vi = vi, Psemax = p, Psemin = 0.999, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        init.vals = NULL)))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lnOR_COPAS2000) = paste0("p = ", p_sa)
lnOR_COPAS2000

#' Copas-N-type selection model (2000)
lnOR_COPAS1999 = vapply(
  p_sa, 
  function(p) {
    mod = suppressWarnings(
      COPAS1999(
        yi = yi, ni = ni, Pnmax = 0.999, Pnmin = p, 
        parset = list(
          mu.bound = 10,
          tau.bound = 5,
          estimate.rho = TRUE, 
          eps = 1e-3,
          init.vals = NULL))
      )
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lnOR_COPAS1999) = paste0("p = ", p_sa)
lnOR_COPAS1999

#' Expected number of missing studies
M_propos = sapply(p_sa, function(p) {
  
  P_max = 0.999
  P_min = p
  ni = n1+n0
  
  n_min = min(ni) 
  n_max = max(ni)
  
  a1 = (qnorm(P_max)-qnorm(P_min))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(P_max)-a1*sqrt(n_max)
  sum((1 - pnorm(a0+a1*sqrt(ni)))/pnorm(a0+a1*sqrt(ni))) 
  
})


M_copas = sapply(p_sa, function(p) {
  
  P_max = p
  P_min = 0.999
  ni = n1+n0
  
  se_min_inv = 1/sqrt(min(vi)) 
  se_max_inv = 1/sqrt(max(vi))
  
  gamma1 = (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 = qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi)))/pnorm(gamma0+gamma1/sqrt(vi)))
  
})


#' Table 1: results of sensitivity analysis
tab1 = data.frame(
  HN = t(lnOR_COPAS_HNGLMM[c(1,3,4),]),
  BN = t(lnOR_COPAS_BNGLMM[c(1,3,4),]),
  CH = t(lnOR_COPAS2000[c(1,3,4),]),  # copas-heckman
  CN = t(lnOR_COPAS1999[c(1,3,4),]),  # copas-n
  M.c = round(M_copas), 
  M.p = round(M_propos)
)

# tab1_p1 = c(
#   HN = lnOR_hn$mu[1], HN.mu.lb = unname(lnOR_hn_lb), HN.mu.ub = unname(lnOR_hn_ub),
#   BN = lnOR_bn$mu[1], BN.mu.lb = unname(lnOR_bn_lb), BN.mu.ub = unname(lnOR_bn_ub),
#   CH = lnOR_nn$mu[1], NN.mu.lb = unname(lnOR_nn_lb), NN.mu.ub = unname(lnOR_nn_ub),
#   CN = lnOR_nn$mu[1], NN.mu.lb = unname(lnOR_nn_lb), NN.mu.ub = unname(lnOR_nn_ub),
#   M.c = 0, M.p = 0
# )

# tab1_all = rbind("p = 1" = tab1_p1, tab1)
tab1_all = tab1
tab1_all$pnmin = c(p_sa)
tab1_all$pnmax = c(rep(0.999, 10))

#' SAVE RESULTS1
#' 
save(lnOR_COPAS_HNGLMM, lnOR_COPAS_BNGLMM, lnOR_COPAS2000, lnOR_COPAS1999,
     M_propos, M_copas, tab1_all,
     lnOR_hn, lnOR_bn, lnOR_nn,
     file = "example-bias.RData")
