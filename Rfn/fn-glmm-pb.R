#' Proposed methods for PB
#' HN-GLMM based sensitivity analysis for lnOR -------
#'
#'@param y0 event in the control group
#'@param y1 event in the experimental group
#'@param n0 subjects in the control group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                estimate.rho: whether to estimate rho    
#'                rho.fix: if not estimate rho, then set a value
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#'   
COPAS_HNGLMM = function(
    y0, y1, n0, n1, 
    Pnmax = 0.99, Pnmin = 0.5,
    n_min,n_max,
    parset = list(
      mu.bound = 10,
      tau.bound = 5,
      estimate.rho = TRUE, 
      rho.fix = -0.9,
      eps = 1e-3,
      integ.limit = 10, 
      cub.tol = 1e-10,
      init.vals = NULL)
    ){
  
  ni = n1+n0
  yi = y1+y0
  
  # n_min = min(ni) 
  # n_max = max(ni)
  
  a1 = (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(Pnmax)-a1*sqrt(n_max)
  
  # if(is.null(parset[["init.vals"]])) {
  #   cyi = metafor::escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0)[,1]
  #   parset[["init.vals"]] = 0.5*c(min(mean(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
  #                             min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
  #   }
  ## estimate rho
  if (parset[["estimate.rho"]]) {

    llk.est.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      rho  = par[3]
      
      f = function(thetai) {
        
        sapply(1:length(yi), 
               function(i) pnorm((a0+a1*sqrt(ni[i])+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(ni[i]))*
                 MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(thetai))*
                # BiasedUrn::dFNCHypergeo(x = y1[i], m1=n1[i], m2=n0[i], n=yi[i], odds=exp(thetai))* 
                 dnorm(thetai, mean = mu, sd = tau)
        )
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(yi), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    ## optimization
    # rho.bound = 0.9999
    # rho.init = -0.1
    
    init.vals.rho = c(parset[["init.vals"]])
    optim.res = try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]], -0.99),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]], 0.99)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      rho  = optim.res$par[3]
      
      hes = numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,9),3,3))
      mu.se  = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      rho.se = sqrt(var.matrix[3,3])
      
      
    } else mu = mu.se = tau2 = tau = tau.se = rho = rho.se = NA
    
    res = list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals.rho,
                var.mat = var.matrix)
    
  } else { ## not estimate rho
    
    rho  = parset[["rho.fix"]]
    llk.fix.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      
      f = function(thetai) {
        
        sapply(1:length(yi), 
               function(i) pnorm((a0+a1*sqrt(ni[i])+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(ni[i]))*
                 MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(thetai))*
                 # BiasedUrn::dFNCHypergeo(x = y1[i], m1=n1[i], m2=n0[i], n=yi[i], odds=exp(thetai)) * 
                 dnorm(thetai, mean = mu, sd = tau)
        )
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(yi), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    optim.res = try(
      nlminb(parset[["init.vals"]], llk.fix.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]]),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      
      hes = numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
      mu.se   = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      
      
    } else mu = mu.se = tau2 = tau = tau.se = NA
    
    res = list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho.fix = parset[["rho.fix"]], rho.se=NA),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = parset[["init.vals"]],
                var.mat = var.matrix)
    
  }
  
  return(res)
  
}

#'
#' BN-GLMM based sensitivity analysis for lnOR -------
#'@param y0 event in the control group
#'@param y1 event in the experimental group
#'@param n0 subjects in the control group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                estimate.rho: whether to estimate rho    
#'                rho.fix: if not estimate rho, then set a value
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#'  
COPAS_BNGLMM = function(
    y0, y1, n0, n1, 
    Pnmax = 0.99, Pnmin = 0.5,
    n_min,n_max,
    parset = list(
      mu.bound = 10,
      tau.bound = 5,
      estimate.rho = TRUE, 
      rho.fix = -0.9,
      eps = 1e-3,
      integ.limit = 10, 
      cub.tol = 1e-10,
      init.vals = NULL)
    ){
  
  ni = n1+n0
  yi = y1+y0

  n_min = min(ni) 
  n_max = max(ni)
  
  a1 = (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(Pnmax)-a1*sqrt(n_max)
  
  if(is.null(parset[["init.vals"]])) {
    cyi = metafor::escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0)[,1]
    parset[["init.vals"]] = 0.5*c(min(mean(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
                              min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
    }
  ## estimate rho
  if (parset[["estimate.rho"]]) {
    
    llk.est.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      rho  = par[3]
      
      f = function(thetai) {
        
        pnorm((a0+a1*sqrt(ni)+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(ni)) * 
          dbinom(y1, yi, prob = plogis(log(n1 / n0) + thetai)) * 
          dnorm(thetai, mean = mu, sd = tau)
        
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(yi), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    ## optimization
    # rho.bound = 0.9999
    # rho.init = -0.1
    
    init.vals.rho = c(parset[["init.vals"]])
    optim.res = try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]], -0.99),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]], 0.99)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      rho  = optim.res$par[3]
      
      hes = numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,9),3,3))
      mu.se  = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      rho.se = sqrt(var.matrix[3,3])
      
      
    } else mu = mu.se = tau2 = tau = tau.se = rho = rho.se = NA
    
    res = list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals.rho,
                var.mat = var.matrix)
    
  } else { ## not estimate rho
    
    rho  = parset[["rho.fix"]]
    llk.fix.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      
      
      f = function(thetai) {
        
        pnorm((a0+a1*sqrt(ni)+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(ni)) * 
          dbinom(y1, yi, prob = plogis(log(n1 / n0) + thetai)) * 
          dnorm(thetai, mean = mu, sd = tau)
        
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(yi), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    ## optimization
    optim.res = try(
      nlminb(parset[["init.vals"]], llk.fix.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]]),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      
      hes = numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
      mu.se   = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])

    } else mu = mu.se = tau2 = tau = tau.se = NA
    
    res = list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho.fix = parset[["rho.fix"]], rho.se=NA),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = parset[["init.vals"]],
                var.mat = var.matrix)
    
  }
  
  return(res)
  
}

#'
#' BN-GLMM based sensitivity analysis for logit-proportion -------
#'@param y0 event in the control group
#'@param y1 event in the experimental group
#'@param n0 subjects in the control group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                estimate.rho: whether to estimate rho    
#'                rho.fix: if not estimate rho, then set a value
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#'  
COPAS_BNGLMM_prop = function(
    y1, n1, 
    Pnmax = 0.99, Pnmin = 0.5,
    n_min,n_max,
    parset = list(
      mu.bound = 10,
      tau.bound = 5,
      estimate.rho = TRUE, 
      rho.fix = -0.9,
      eps = 1e-3,
      integ.limit = 10, 
      cub.tol = 1e-10,
      init.vals = NULL)
    ){
  
  n_min = min(n1) 
  n_max = max(n1)
  
  a1 = (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(Pnmax)-a1*sqrt(n_max)

  # if(is.null(parset[["init.vals"]])) {
  #   cyi = metafor::escalc(measure="PLN", xi=y1, mi=n1-y1)[,1]
  #   parset[["init.vals"]] = 0.5*c(min(median(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
  #                             min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
  #   }

  ## estimate rho
  if (parset[["estimate.rho"]]) {
    
    llk.est.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      rho  = par[3]
      
      f = function(thetai) {
        
        pnorm((a0+a1*sqrt(n1)+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(n1)) * 
          dbinom(y1, n1, prob = plogis(thetai)) * 
          dnorm(thetai, mean = mu, sd = tau)
        
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(y1), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }

    ## optimization
    # rho.bound = 0.9999
    # rho.init = -0.1
    
    init.vals.rho = c(parset[["init.vals"]])
    optim.res = try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]], -0.99),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]], 0.99)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      rho  = optim.res$par[3]
      
      hes = numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,9),3,3))
      mu.se  = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      rho.se = sqrt(var.matrix[3,3])
      
      
    } else mu = mu.se = tau2 = tau = tau.se = rho = rho.se = NA
    
    res = list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals.rho,
                var.mat = var.matrix)
    
  } else {  ## not estimate rho
    
    
    llk.fix.rho = function(par) {
      
      mu   = par[1]
      tau  = par[2]
      tau2 = tau^2
      rho  = parset[["rho.fix"]]
      
      f = function(thetai) {
        
        pnorm((a0+a1*sqrt(n1)+rho*(thetai-mu)/tau)/suppressWarnings(sqrt(1-rho^2)))/pnorm(a0+a1*sqrt(n1)) * 
          dbinom(y1, n1, prob = plogis(thetai)) * 
          dnorm(thetai, mean = mu, sd = tau)
        
      }
      
      prob.prior = cubature::hcubature(f, 
        lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
        fDim = length(yi), tol = parset[["cub.tol"]])$integral
      l = sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }

    ## optimization
    optim.res = try(
      nlminb(parset[["init.vals"]], llk.fix.rho,
             lower = c(-parset[["mu.bound"]], parset[["eps"]]),
             upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   = optim.res$par[1]
      tau  = optim.res$par[2]
      tau2 = tau^2
      
      hes = numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] = sqrt(parset[["eps"]])
      var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
      mu.se   = sqrt(var.matrix[1,1])
      tau.se = sqrt(var.matrix[2,2])
      
    } else mu = mu.se = tau2 = tau = tau.se = NA
    
    res = list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho.fix = parset[["rho.fix"]], rho.se=NA),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = var.matrix)
    
  }
  
  return(res)
  
}

