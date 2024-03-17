##
## COPAS-HECKMAN-TYPE SELECTION MODEL
##
## 

lambda <- function(x)
  dnorm(x) / pnorm(x)

## LIKELIHOOD ESTIMATING rho, FROM metasens::copas()
copas.loglik.without.beta.est.rho <- function(
    x, gamma,
    TE, seTE) {
  
  
  mu  <- x[1]
  rho <- x[3]
  tau <- x[2]
  ##
  ## TE   <=> estimated treatment effect
  ## seTE <=> standard error from trials, conditional on publication
  
  
  ## Copas, Shi (2000), Biostatistics, p. 250:
  ##
  u <- gamma[1] + gamma[2] / seTE
  ##
  sigma <- sqrt(seTE^2 / (1 - rho^2 * lambda(u) * (u + lambda(u))))
  ##
  s2t2 <- sigma^2 + tau^2
  ##
  rho.tilde <- rho * sigma / sqrt(s2t2)
  ##
  v <- (u + rho.tilde * (TE - mu) / (sqrt(s2t2))) /
    sqrt(1 - rho.tilde^2)
  ##
  ## Avoid numerical problems by replacing 0's in pnorm(v):
  ## qnorm(1e-320) = -38.26913
  ## this is towards the smallest value for log
  ##
  v[v < -37] <- -37
  ##
  ## Take minus log-likelihood and minimise it;
  ## leave out log(pnorm(u)) as this is a constant
  ##
  ell <- -(-0.5 * log(s2t2) - (TE - mu)^2 / (2 * s2t2) + log(pnorm(v)))
  
  res <- sum(ell)
  ##
  res
}

## LIKELIHOOD NOT ESTIMATING rho
copas.loglik.without.beta.fix.rho <- function(
    x, gamma, rho,
    TE, seTE) {
  
  
  mu  <- x[1]
  tau <- x[2]
  ##
  ## TE   <=> estimated treatment effect
  ## seTE <=> standard error from trials, conditional on publication
  
  
  ## Copas, Shi (2000), Biostatistics, p. 250:
  ##
  u <- gamma[1] + gamma[2] / seTE
  ##
  sigma <- sqrt(seTE^2 / (1 - rho^2 * lambda(u) * (u + lambda(u))))
  ##
  s2t2 <- sigma^2 + tau^2
  ##
  rho.tilde <- rho * sigma / sqrt(s2t2)
  ##
  v <- (u + rho.tilde * (TE - mu) / (sqrt(s2t2))) /
    sqrt(1 - rho.tilde^2)
  ##
  ## Avoid numerical problems by replacing 0's in pnorm(v):
  ## qnorm(1e-320) = -38.26913
  ## this is towards the smallest value for log
  ##
  v[v < -37] <- -37
  ##
  ## Take minus log-likelihood and minimise it;
  ## leave out log(pnorm(u)) as this is a constant
  ##
  ell <- -(-0.5 * log(s2t2) - (TE - mu)^2 / (2 * s2t2) + log(pnorm(v)))
  
  res <- sum(ell)
  ##
  res
}


##
## COPAS-HECKMAN-TYPE SELECTION MODEL
##
copas_NNLMM <- function(
    yi, vi, Psemax = 0.99, Psemin = 0.5,
    ## IF rho.sa = NULL, THEN ESTIMATE rho, AND rho.upper, rho.init WORK
    rho.sa = NULL, ## 
    rho.upper = 0.9999, rho.init = -0.1,
    mu.upper = Inf,
    tau.upper = 3,
    eps = .Machine$double.eps,
    init.vals = c(lnOR_nn$mu[1], lnOR_nn$tau[1])){
  
  se <- sqrt(vi)
  # se_min <- min(se) 
  # se_max <- max(se)
  # 
  # gamma1 <- (qnorm(Psemax)-qnorm(Psemin))/(1/se_maxn-1/se_min)
  # gamma0 <- qnorm(Psemax)-gamma1/se_max
  # 
  se_min_inv <- 1/min(se)
  se_max_inv <- 1/max(se)
  
  gamma1 <- (qnorm(Psemax)-qnorm(Psemin))/(se_max_inv-se_min_inv)
  gamma0 <- qnorm(Psemax)-gamma1*se_max_inv
  
  gamma <- c(gamma0, gamma1)
  
  if (is.null(rho.sa)) {
    ## estimate rho
    
    
    if(is.null(init.vals)) init.vals <- c(mu = mean(yi, na.rm = TRUE), tau = mean(vi, na.rm = TRUE))
    init.vals.rho <- c(init.vals, rho = rho.init)
    
    llk.est.rho <- function(x) copas.loglik.without.beta.est.rho(x, gamma = gamma, TE = yi, seTE = sqrt(vi))
    optim.res <- try(
      nlminb(init.vals.rho, llk.est.rho, ## mu tau rho
             lower = c(-mu.upper, eps, -rho.upper),
             upper = c( mu.upper, tau.upper,  rho.upper)), 
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
                gamma = gamma,
                opt = optim.res,
                init.vals = init.vals.rho)
    
  } else {
    ## fix value of rho
    
    if(is.null(init.vals)) init.vals <- c(mu = mean(yi, na.rm = TRUE), tau = mean(vi, na.rm = TRUE))
    llk.fix.rho <- function(x) copas.loglik.without.beta.fix.rho(x, gamma = gamma, rho = rho.sa, TE = yi, seTE = sqrt(vi))
    optim.res <- try(
      nlminb(init.vals, llk.fix.rho, ## mu tau
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
                gamma = gamma,
                opt = optim.res,
                init.vals = init.vals)
    
  }
  
  return(res)
  
}

# copas_NNLMM(yi = yi,vi = vi, Psemax = 0.99, Psemin = 0.1, rho.init = -0.1)$opt$par
# copas_NNLMM(yi = yi,vi = vi, Psemax = 0.99, Psemin = 0.1, rho.sa = 0.9999)$opt$par




