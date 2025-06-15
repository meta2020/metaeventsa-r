#' Copa-Heckman-type selection model for adjusting PB
#' [1] COPAS J, SHI J. Meta-analysis, funnel plots and sensitivity analysis. 
#' Biostatistics, 2000, 1(3): 247-262. 

.lambda = function(x) dnorm(x) / pnorm(x)
#' COPAS-HECKMAN-TYPE SELECTION MODEL
#' 
#'@param yi outcomes
#'@param vi variance of the outcomes (SE^2)
#'@param Psemax probabilities of publishing a study with maximum SE 
#'@param Psemin probabilities of publishing a study with minimum SE 
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                estimate.rho: whether to estimate rho    
#'                rho.fix: if not estimate rho, then set a value
#'                eps: smallest value used
#'                init.vals: initial values of mu and tau
#'
COPAS2000 = function(
    yi, vi, 
    Psemax = 0.99, Psemin = 0.5,
    parset = list(
      mu.bound = 10,
      tau.bound = 5,
      estimate.rho = TRUE, 
      rho.fix = -0.9,
      eps = 1e-3,
      init.vals = NULL)
    ){
  
  si = sqrt(vi)
  # se_min = min(se) 
  # se_max = max(se)
  # 
  # gamma1 = (qnorm(Psemax)-qnorm(Psemin))/(1/se_maxn-1/se_min)
  # gamma0 = qnorm(Psemax)-gamma1/se_max
  # 
  se_min_inv = 1/min(si)
  se_max_inv = 1/max(si)
  
  gamma1 = (qnorm(Psemax)-qnorm(Psemin))/(se_max_inv-se_min_inv) ## gamma1 tends to be negative
  gamma0 = qnorm(Psemax)-gamma1*se_max_inv
  
  gamma = c(gamma0, gamma1)

  # if(is.null(parset[["init.vals"]])) parset[["init.vals"]] = c(median(yi, na.rm = TRUE), sd(yi, na.rm = TRUE))

  ## estimate rho
  if (parset[["estimate.rho"]]) {

    copas2000.llk.est.rho = function(par) {
  
      mu  = par[1]
      tau = par[2]
      tau2= tau^2
      rho = par[3]
      
      ## Copas, Shi (2000), Biostatistics, p. 250:
      ##
      u = gamma[1] + gamma[2] / si
      ##
      sigma = sqrt(si^2 / (1 - rho^2 * .lambda(u) * (u + .lambda(u))))
      s2t2 = sigma^2 + tau^2
      rho.tilde = rho * sigma / sqrt(s2t2)
      v = (u + rho.tilde * (yi - mu) / (sqrt(s2t2))) / suppressWarnings(sqrt(1 - rho.tilde^2))
      ##
      ## Avoid numerical problems by replacing 0's in pnorm(v):
      ## qnorm(1e-320) = -38.26913
      ## this is towards the smallest value for log
      ##
      v[v < -37] = -37
      ##
      ## Take minus log-likelihood and minimise it;
      ## leave out log(pnorm(u)) as this is a constant
      ##
      ell = -(-0.5 * log(s2t2) - (yi - mu)^2 / (2 * s2t2) + log(pnorm(v)))
      
      res = sum(ell)
      ##
      return(res)
    }
    ## optimization
    init.vals.rho = c(parset[["init.vals"]])
    optim.res = try(
      nlminb(init.vals.rho, copas2000.llk.est.rho,
               lower = c(-parset[["mu.bound"]], parset[["eps"]], -0.99),
               upper = c( parset[["mu.bound"]], parset[["tau.bound"]], 0.99)), 
        silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      rho  = optim.res$par[3]
      
      hes = numDeriv::hessian(copas2000.llk.est.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,9),3,3))
      mu.se  = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      rho.se = sqrt(var.matrix[3,3])
      
    } else mu = mu.se = tau2 = tau = tau.se = rho = rho.se = NA
    
    res = list(mu  = c(mu = mu, mu.se = mu.se),
               tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
               rho = c(rho = rho, rho.se = rho.se),
               gamma = gamma,
               opt = optim.res,
               init.vals = init.vals.rho,
               var.mat = var.matrix)
    
  } else {  ## not estimate rho
    
  ## likelihood function 
  copas2000.llk.fix.rho = function(par){
  
    mu  = par[1]
    tau = par[2]
    tau2= tau^2
    rho = rho.fix

    u = gamma[1] + gamma[2] / seTE
    sigma = sqrt(seTE^2 / (1 - rho^2 * .lambda(u) * (u + .lambda(u))))
    s2t2 = sigma^2 + tau^2
    rho.tilde = rho * sigma / sqrt(s2t2)
    ##
    v = (u + rho.tilde * (TE - mu) / (sqrt(s2t2))) /
      suppressWarnings(sqrt(1 - rho.tilde^2))
    ##
    ## Avoid numerical problems by replacing 0's in pnorm(v):
    ## qnorm(1e-320) = -38.26913
    ## this is towards the smallest value for log
    ##
    v[v < -37] = -37
    ##
    ## Take minus log-likelihood and minimise it;
    ## leave out log(pnorm(u)) as this is a constant
    ##
    ell = -(-0.5 * log(s2t2) - (TE - mu)^2 / (2 * s2t2) + log(pnorm(v)) - log(pnorm(u)))
    
    res = sum(ell)
    ##
    return(res)
  }

    ## optimization
    optim.res = try(
      nlminb(parset[["init.vals"]], copas2000.llk.fix.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]]),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      
      hes = numDeriv::hessian(copas2000.llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(eps)
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
      mu.se   = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      
      
    } else mu = mu.se = tau2 = tau = tau.se = NA
    
    res = list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = rho.fix,
                gamma = gamma,
                opt = optim.res,
                init.vals = init.vals,
                var.mat = var.matrix)
    
  }
  
  return(res)
  
}





