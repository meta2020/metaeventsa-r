#' Original models without considering PB
#' NN model for lnOR/logit-proportion -------
#'
#'@param yi estimated lnOR
#'@param vi SE of the lnOR
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                eps: smallest value used
#'                init.vals: initial values of mu and tau
#'               
NN_LMM = function(
    yi, vi, 
    parset=list(
      mu.bound = 10, 
      tau.bound = 5,
      eps = 0.001,
      init.vals = NULL
      )
    ){
  
  ## likelihood function
  llk.fn = function(par) {
    
    mu   = par[1]
    tau  = par[2]
    tau2 = tau^2
   
    llk = dnorm(yi, mean = mu, sd = sqrt(vi + tau2), log = T)
    l = sum(llk, na.rm = TRUE)
    
    return(-l)
  }
  
  ## optimization
  if(is.null(parset[["init.vals"]])) {
    cyi = yi
    parset[["init.vals"]] = 0.5*c(min(mean(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
                              min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
  }
  optim.res = tryCatch(
    nlminb(parset[["init.vals"]], llk.fn,
    lower = c(-parset[["mu.bound"]], parset[["eps"]]),
    upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
    error = function(e) NULL)
  
  ## summarize results
  if(!is.null(optim.res)) {
    
    mu   = optim.res$par[1]
    tau  = optim.res$par[2]
    tau2 = tau^2
    
    hes = numDeriv::hessian(llk.fn, optim.res$par)
    hes[is.nan(hes)] = sqrt(parset[["eps"]])
    var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
    mu.se   = sqrt(var.matrix[1,1])
    tau.se = sqrt(var.matrix[2,2])
    
  } else mu = mu.se = tau2 = tau = tau.se = NA
  
  res = list(mu = c(mu = mu, mu.se = mu.se),
             tau = c(tau = tau, tau.se = tau.se, tau2 = tau2) ,
             opt = optim.res,
             init.vals = parset[["init.vals"]],
             var.mat = var.matrix)
  
  return(res)
  
  
}

#' HN-GLMM for lnOR -------
#'
#'@param y0 event in the control group
#'@param y1 event in the experimental group
#'@param n0 subjects in the control group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#'   
HN_GLMM = function(
    y0, y1, n0, n1,
    parset = list(
      mu.bound = 10,
      tau.bound = 3,
      eps = 0.001,
      integ.limit = 10, 
      cub.tol = 1e-5,
      init.vals = NULL)
    ){
  
  yi = y0+y1
  
  ## likelihood function
  llk.fn = function(par) {
    
    mu  = par[1]
    tau = par[2]
    
    f = function(thetai) {  
      sapply(1:length(yi), 
        function(i) MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(thetai)) *
        # function(i) BiasedUrn::dFNCHypergeo(x = y1[i], m1=n1[i], m2=n0[i], n=yi[i], odds=exp(thetai)) * 
          dnorm(thetai, mean = mu, sd = tau))
      }

    prob.prior = cubature::hcubature(f, 
      lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
      fDim = length(yi), tol = parset[["cub.tol"]])$integral
    l = sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }

  ## optimization
  if(is.null(parset[["init.vals"]])) {
    cyi = metafor::escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0)[,1]
    parset[["init.vals"]] = 0.5*c(min(mean(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
                              min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
    }
  optim.res = try(
    nlminb(parset[["init.vals"]], llk.fn, 
           lower = c(-parset[["mu.bound"]], parset[["eps"]]),
           upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
    silent = TRUE)
  
  ## summarize results
  if(!inherits(optim.res, "try-error")) {
    
    mu   = optim.res$par[1]
    tau  = optim.res$par[2]
    tau2 = tau^2
    
    hes = numDeriv::hessian(llk.fn, optim.res$par)
    hes[is.nan(hes)] = sqrt(parset[["eps"]])
    var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
    mu.se  = sqrt(var.matrix[1,1])
    tau.se = sqrt(var.matrix[2,2])
    
  } else mu = mu.se = tau = tau2 = tau.se = NA
  
  res = list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = parset[["init.vals"]],
             var.mat = var.matrix)
  
  return(res)
  
  
}

#' BN-GLMM for lnOR -------
#'
#'@param y0 event in the control group
#'@param y1 event in the experimental group
#'@param n0 subjects in the control group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#'     
BN_GLMM = function(
    y0, y1, n0, n1,
    parset = list(
      mu.bound = 10,
      tau.bound = 3,
      eps = 0.001,
      integ.limit = 10, 
      cub.tol = 1e-5,
      init.vals = NULL)
    ){
  
  yi = y0+y1
  
  ## likelihood function
  llk.fn = function(par) {
    
    mu = par[1]
    tau = par[2]
    
    f = function(thetai) dbinom(
      y1, yi, prob = plogis(log(n1 / n0) + thetai)) * dnorm(thetai, mean = mu, sd = tau)
    
    prob.prior = cubature::hcubature(f, 
      lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], fDim = length(y1), tol = parset[["cub.tol"]])$integral
    l = sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }

  ## optimization
  if(is.null(parset[["init.vals"]])) {
    cyi = metafor::escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0)[,1]
    parset[["init.vals"]] = 0.5*c(min(mean(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
                              min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
    }
  optim.res = try(
    nlminb(parset[["init.vals"]], llk.fn,
           lower = c(-parset[["mu.bound"]], parset[["eps"]]),
           upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
    silent = TRUE)

  ## summarize results
  if(!inherits(optim.res, "try-error")) {
    
    mu   = optim.res$par[1]
    tau  = optim.res$par[2]
    tau2 = tau^2
    
    hes = numDeriv::hessian(llk.fn, optim.res$par)
    hes[is.nan(hes)] = sqrt(parset[["eps"]])
    var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
    mu.se  = sqrt(var.matrix[1,1])
    tau.se = sqrt(var.matrix[2,2])
  
    } else mu = mu.se = tau = tau2 = tau.se = NA
 
  res = list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = parset[["init.vals"]],
             var.mat = var.matrix)
  
  return(res)
  
}

#' BN-GLMM for logit-proportion -------
#'
#'@param y1 event in the experimental group
#'@param n1 subjects in the experimental group
#'@param parset additional parameters setting: 
#'                mu.bound: upper and lower(negative value) bounds for mu
#'                tau.bound: upper bound for tau
#'                eps: smallest value used
#'                integ.limit: upper and lower limit in double integral
#'                cub.tol: tolerance in double integral
#'                init.vals: initial values of mu and tau
#' 
BN_GLMM_prop = function(
    y1, n1,
    parset = list(
      mu.bound = 10,
      tau.bound = 3,
      eps = 0.001,
      integ.limit = 10, 
      cub.tol = 1e-5,
      init.vals = NULL)
    ){
  ## likelihood function
  llk.fn = function(par) {
    
    mu = par[1]
    tau = par[2]
    
    f = function(thetai) dbinom(
      y1, n1, prob = plogis(thetai)) * dnorm(thetai, mean = mu, sd = tau)
    
    prob.prior =  cubature::hcubature(f, 
      lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
      fDim = length(y1), tol = parset[["cub.tol"]])$integral
    l = sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }

  ## optimization
    if(is.null(parset[["init.vals"]])) {
    cyi = metafor::escalc(measure="PLN", xi=y1, mi=n1-y1)[,1]
    parset[["init.vals"]] = 0.5*c(min(median(cyi, na.rm = TRUE), parset[["mu.bound"]]), 
                              min(sd(cyi, na.rm = TRUE), parset[["tau.bound"]]))
    }
  optim.res = try(
    nlminb(parset[["init.vals"]], llk.fn,
           lower = c(-parset[["mu.bound"]], parset[["eps"]]),
           upper = c( parset[["mu.bound"]], parset[["tau.bound"]])), 
    silent = TRUE)
  
  if(!inherits(optim.res, "try-error")) {
    
    mu   = optim.res$par[1]
    tau  = optim.res$par[2]
    tau2 = tau^2
    
    hes = numDeriv::hessian(llk.fn, optim.res$par)
    hes[is.nan(hes)] = sqrt(parset[["eps"]])
    var.matrix = tryCatch(solve(hes), error=function(e) matrix(rep(NA,4),2,2))
    mu.se  = sqrt(var.matrix[1,1])
    tau.se = sqrt(var.matrix[2,2])
    
  } else mu = mu.se = tau = tau2 = tau.se = NA
  
  res = list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = parset[["init.vals"]],
             var.mat = var.matrix)
  
  return(res)
  
}


