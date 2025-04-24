##
## HN-GLMM BASED SENSITIVITY ANALYSIS FOR lnOR ---------------------------------
##

copas_HNGLMM <- function(
    y0, y1, n0, n1, Pnmax = 0.99, Pnmin = 0.5,
    ## IF rho.sa = NULL, THEN ESTIMATE rho, AND rho.upper, rho.init WORK
    rho.sa = NULL, ## 
    rho.upper = 0.9999, rho.init = -0.1,
    mu.upper = 3,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = 10, cub.tol = 1e-10,
    init.vals = c(lnOR_bn$mu[1], lnOR_bn$tau[1])){
  
  ni <- n1+n0
  yi <- y1+y0
  
  n_min <- min(ni) 
  n_max <- max(ni)
  
  a1 <- (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 <- qnorm(Pnmax)-a1*sqrt(n_max)
  
  if (is.null(rho.sa)) {
    ## estimate rho
    
    llk.est.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- par[3]
      
      f <- function(theta) {
        
        sapply(1:length(yi), 
               function(i) pnorm((a0+a1*sqrt(ni[i])+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(ni[i]))*
                 MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(theta))* 
                 dnorm(theta, mean = mu, sd = tau)
        )
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(yi), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    init.vals.rho <- c(init.vals, rho = rho.init)
    
    optim.res <- try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-mu.upper, eps, -rho.upper),
             upper = c( mu.upper, tau.upper, rho.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      rho  <- optim.res$par[3]
      
      hes <- numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se  <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      rho.se <- sqrt(var.matrix[3,3])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- rho <- rho.se <- NA
    
    res <- list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  } else {
    ## fix value of rho
    
    llk.fix.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- rho.sa
      
      f <- function(theta) {
        
        sapply(1:length(yi), 
               function(i) pnorm((a0+a1*sqrt(ni[i])+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(ni[i]))*
                 MCMCpack::dnoncenhypergeom(x = y1[i], n1[i], n0[i], yi[i], exp(theta))* 
                 dnorm(theta, mean = mu, sd = tau)
        )
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(yi), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    if(is.null(init.vals)) init.vals <- c(mu = mean(yi, na.rm = TRUE), tau = mean(vi, na.rm = TRUE))
    
    optim.res <- try(
      nlminb(init.vals, llk.fix.rho,
             lower = c(-mu.upper, eps),
             upper = c( mu.upper, tau.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      
      hes <- numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se   <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- NA
    
    res <- list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = rho.sa,
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  }
  
  return(res)
  
}

# copas_HNGLMM(y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.9, Pnmin = 0.9, rho.init = 0.1, init.vals = c(lnOR_hn$mu[1], lnOR_hn$tau[1]))
# copas_HNGLMM(y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.99, Pnmin = 0.5, rho.sa = -0.99)




##
## BN-GLMM BASED SENSITIVITY ANALYSIS FOR lnOR ---------------------------------
##

copas_BNGLMM <- function(
    y0, y1, n0, n1, Pnmax = 0.99, Pnmin = 0.5,
    ## IF rho.sa = NULL, THEN ESTIMATE rho, AND rho.upper, rho.init WORK
    rho.sa = NULL, ## 
    rho.upper = 0.9999, rho.init = -0.1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = 10, cub.tol = 1e-10,
    init.vals = c(lnOR_bn$mu[1], lnOR_bn$tau[1])){
  
  ni <- n1+n0
  yi <- y1+y0

  n_min <- min(ni) 
  n_max <- max(ni)
  
  a1 <- (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 <- qnorm(Pnmax)-a1*sqrt(n_max)
  
  if (is.null(rho.sa)) {
    ## estimate rho
    
    llk.est.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- par[3]
      
      f <- function(theta) {
        
        pnorm((a0+a1*sqrt(ni)+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(ni)) * 
          dbinom(y1, yi, prob = plogis(log(n1 / n0) + theta)) * 
          dnorm(theta, mean = mu, sd = tau)
        
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(yi), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    init.vals.rho <- c(init.vals, rho = rho.init)
    
    optim.res <- try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-mu.upper, eps, -rho.upper),
             upper = c( mu.upper, tau.upper, rho.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      rho  <- optim.res$par[3]
      
      hes <- numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se  <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      rho.se <- sqrt(var.matrix[3,3])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- rho <- rho.se <- NA
    
    res <- list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  } else {
    ## fix value of rho
    
    llk.fix.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- rho.sa
      
      f <- function(theta) {
        
        pnorm((a0+a1*sqrt(ni)+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(ni)) * 
          dbinom(y1, yi, prob = plogis(log(n1 / n0) + theta)) * 
          dnorm(theta, mean = mu, sd = tau)
        
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(yi), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    if(is.null(init.vals)) init.vals <- c(mu = mean(yi, na.rm = TRUE), tau = mean(vi, na.rm = TRUE))
    
    optim.res <- try(
      nlminb(init.vals, llk.fix.rho,
             lower = c(-mu.upper, eps),
             upper = c( mu.upper, tau.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      
      hes <- numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se   <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- NA
    
    res <- list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = rho.sa,
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  }
  
  return(res)
  
}

# copas_BNGLMM(y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.99, Pnmin = 0.1, rho.init = -0.1)$opt$par
# copas_BNGLMM(y0=y0, y1=y1, n0=n0, n1=n1, Pnmax = 0.99, Pnmin = 0.1, rho.sa = -0.2)$opt$par

##
## BN-GLMM BASED SENSITIVITY ANALYSIS FOR logit-proportion ---------------------
##

copas_BNGLMM_prop <- function(
    y1, n1, Pnmax = 0.99, Pnmin = 0.5,
    ## IF rho.sa = NULL, THEN ESTIMATE rho, AND rho.upper, rho.init WORK
    rho.sa = NULL, ## 
    rho.upper = 0.9999, rho.init = -0.1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    integ.limit = 10, cub.tol = 1e-10,
    init.vals = c(lnOR_bn$mu[1], lnOR_bn$tau[1])){
  
  n_min <- min(n1) 
  n_max <- max(n1)
  
  a1 <- (qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0 <- qnorm(Pnmax)-a1*sqrt(n_max)
  
  if (is.null(rho.sa)) {
    ## estimate rho
    
    llk.est.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- par[3]
      
      f <- function(theta) {
        
        pnorm((a0+a1*sqrt(n1)+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(n1)) * 
          dbinom(y1, n1, prob = plogis(theta)) * 
          dnorm(theta, mean = mu, sd = tau)
        
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(y1), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }
    
    init.vals.rho <- c(init.vals, rho = rho.init)
    
    optim.res <- try(
      nlminb(init.vals.rho, llk.est.rho,
             lower = c(-mu.upper, eps, -rho.upper),
             upper = c( mu.upper, tau.upper, rho.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      rho  <- optim.res$par[3]
      
      hes <- numDeriv::hessian(llk.est.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se  <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      rho.se <- sqrt(var.matrix[3,3])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- rho <- rho.se <- NA
    
    res <- list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = c(rho = rho, rho.se = rho.se),
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  } else {
    ## fix value of rho
    
    llk.fix.rho <- function(par) {
      
      mu   <- par[1]
      tau  <- par[2]
      tau2 <- tau^2
      rho  <- rho.sa
      
      f <- function(theta) {
        
        pnorm((a0+a1*sqrt(n1)+rho*(theta-mu)/tau)/sqrt(1-rho^2))/pnorm(a0+a1*sqrt(n1)) * 
          dbinom(y1, n1, prob = plogis(theta)) * 
          dnorm(theta, mean = mu, sd = tau)
        
      }
      
      prob.prior <- cubature::hcubature(f, lowerLimit = -integ.limit, upperLimit = integ.limit, fDim = length(y1), tol = cub.tol)$integral
      
      l <- sum(log(prob.prior), na.rm = TRUE)
      
      return(-l)
    }

    optim.res <- try(
      nlminb(init.vals, llk.fix.rho,
             lower = c(-mu.upper, eps),
             upper = c( mu.upper, tau.upper)), 
      silent = TRUE)
    
    if(!inherits(optim.res, "try-error")) {
      
      mu   <- optim.res$par[1]
      tau  <- optim.res$par[2]
      tau2 <- tau^2
      
      hes <- numDeriv::hessian(llk.fix.rho, optim.res$par)
      hes[is.nan(hes)] <- sqrt(eps)
      var.matrix <- solve(hes)
      mu.se   <- sqrt(var.matrix[1,1])
      tau.se <- sqrt(var.matrix[2,2])
      
      
    } else mu <- mu.se <- tau2 <- tau <- tau.se <- NA
    
    res <- list(mu = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                rho = rho.sa,
                a   = c(a0, a1),
                opt = optim.res,
                init.vals = init.vals,
                var.mat = solve(hes))
    
  }
  
  return(res)
  
}

