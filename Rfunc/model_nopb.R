##
## NN MODEL FOR lnOR/logit-proportion ------------------------------------------
##
NN_LMM <- function(
    yi, vi, 
    init.vals = c(mean(yi, na.rm = TRUE), mean(vi, na.rm = TRUE)), 
    mu.upper = Inf, tau.upper = 3,
    eps = .Machine$double.eps){
  
  
  llk.fn = function(par) {
    
    mu   <- par[1]
    tau  <- par[2]
    tau2 <- tau^2
   
    llk <- dnorm(yi, mean = mu, sd = sqrt(vi + tau2), log = T)
    
    l <- sum(llk, na.rm = TRUE)
    
    return(-l)
  }
  
  optim.res <- try(
    nlminb(init.vals, llk.fn,
    lower = c(-mu.upper, eps),
    upper = c( mu.upper, tau.upper)), silent = TRUE)
  
  if(!inherits(optim.res, "try-error")) {
    
    mu   <- optim.res$par[1]
    tau  <- optim.res$par[2]
    tau2 <- tau^2
    
    hes <- numDeriv::hessian(llk.fn, optim.res$par)
    var.matrix <- solve(hes)
    mu.se   <- sqrt(var.matrix[1,1])
    tau.se <- sqrt(var.matrix[2,2])
    
    
  } else mu <- mu.se <- tau2 <- tau <- tau.se <- NA
  
  res <- list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2) ,
              opt = optim.res,
              init.vals = init.vals)
  
  return(res)
  
  
}

##
## HN-GLMM FOR lnOR ------------------------------------------------------------
##
HN_GLMM <- function(
    y0, y1, n0, n1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = 10, cub.tol = 1e-5,
    init.vals = c(-0.1, 0.1)){
  
  yi <- y0+y1
  
  llk.fn <- function(par) {
    
    mu  <- par[1]
    tau <- par[2]
    
    f <- function(theta) {
      
      sapply(1:length(yi), 
             function(i) MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(theta))* dnorm(theta, mean = mu, sd = tau)
      )
      }
      

    prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(yi), tol = cub.tol)$integral
    
    l <- sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }
  
  optim.res <- try(
    nlminb(init.vals, llk.fn, 
           lower = c(-mu.upper, eps),
           upper = c( mu.upper, tau.upper)), 
    silent = TRUE)
  
  if(!inherits(optim.res, "try-error")) {
    
    mu   <- optim.res$par[1]
    tau  <- optim.res$par[2]
    tau2 <- tau^2
    
    hes <- numDeriv::hessian(llk.fn, optim.res$par)
    var.matrix <- solve(hes)
    mu.se  <- sqrt(var.matrix[1,1])
    tau.se <- sqrt(var.matrix[2,2])
    
  } else mu <- mu.se <- tau <- tau2 <- tau.se <- NA
  
  res <- list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = init.vals)
  
  return(res)
  
  
}



##
## BN-GLMM FOR lnOR ------------------------------------------------------------
##

BN_GLMM <- function(
    y0, y1, n0, n1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = Inf, cub.tol = 1e-5,
    init.vals = c(-0.1, 0.1)){
  
  yi <- y0+y1
  
  llk.fn <- function(par) {
    
    mu <- par[1]
    tau <- par[2]
    
    f <- function(theta) dbinom(y1, yi, prob = plogis(log(n1 / n0) + theta)) * dnorm(theta, mean = mu, sd = tau)
    
    prob.prior <-  cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(y1), tol = cub.tol)$integral
    
    l <- sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }
  
  optim.res <- try(
    nlminb(init.vals, llk.fn,
           lower = c(-mu.upper, eps),
           upper = c( mu.upper, tau.upper)), 
    silent = TRUE)
  
  if(!inherits(optim.res, "try-error")) {
    
    mu   <- optim.res$par[1]
    tau  <- optim.res$par[2]
    tau2 <- tau^2
    
    hes <- numDeriv::hessian(llk.fn, optim.res$par)
    var.matrix <- solve(hes)
    mu.se  <- sqrt(var.matrix[1,1])
    tau.se <- sqrt(var.matrix[2,2])
  
    } else mu <- mu.se <- tau <- tau2 <- tau.se <- NA
 
  res <- list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = init.vals)
  
  return(res)
  
}





##
## BN-GLMM FOR logit-proportion
##

BN_GLMM_prop <- function(
    y1, n1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = Inf, cub.tol = 1e-5,
    init.vals = c(-0.1, 0.1)){
  
  llk.fn <- function(par) {
    
    mu <- par[1]
    tau <- par[2]
    
    f <- function(theta) dbinom(y1, n1, prob = plogis(theta)) * dnorm(theta, mean = mu, sd = tau)
    
    prob.prior <-  cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(y1), tol = cub.tol)$integral
    
    l <- sum(log(prob.prior), na.rm = TRUE)
    
    return(-l)
  }
  
  optim.res <- try(
    nlminb(init.vals, llk.fn,
           lower = c(-mu.upper, eps),
           upper = c( mu.upper, tau.upper)), 
    silent = TRUE)
  
  if(!inherits(optim.res, "try-error")) {
    
    mu   <- optim.res$par[1]
    tau  <- optim.res$par[2]
    tau2 <- tau^2
    
    hes <- numDeriv::hessian(llk.fn, optim.res$par)
    var.matrix <- solve(hes)
    mu.se  <- sqrt(var.matrix[1,1])
    tau.se <- sqrt(var.matrix[2,2])
    
  } else mu <- mu.se <- tau <- tau2 <- tau.se <- NA
  
  res <- list(mu = c(mu = mu, mu.se = mu.se),
              tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
              opt = optim.res,
              init.vals = init.vals)
  
  return(res)
  
}


