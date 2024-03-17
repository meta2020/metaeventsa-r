##
## META-ANALYSIS OF THE logit-proportion
## ESTIMATIONS WIH AND WITHOUT PUBLICATION BIAS
##
## LOAD R FUNCTIONS
file.sources <- list.files("Rfunc/")
sapply(paste0("Rfunc/", file.sources), source)

## DATA ARE FROM Stijnen ET AL, 2010
data <- read.csv("niel-weise21.csv")
data$c1 <- data$n1-data$y1
data$c0 <- data$n0-data$y0

## DATA AFTER CONTINUITY CORRECTION
data_cc <- cc(data)

## ESTIMATION WITHOU PUBLICATION BIAS (Stijnen ET AL, 2010)
### DATA OF EVENTS AND SUBJECTS
y1 <- data$y1
n1 <- data$n1


## BN-GLMM
lnOR_bn <- BN_GLMM_prop(
  y1 = y1, n1 = n1, 
  mu.upper = 5, tau.upper = 3, integ.limit = 10, 
  init.vals = c(-3,0.1))
lnOR_bn_lb <- lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_bn$mu[2]
lnOR_bn_ub <- lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_bn$mu[2]

sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_bn$mu[1], lnOR_bn_lb, lnOR_bn_ub, lnOR_bn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lnOR_bn$tau[1], lnOR_bn$tau[2])

## NN RANDOM-EFFECTS MODEL
### CALCULATE THE lnOR AND SE
yi <- log(data_cc$y1/data_cc$c1)
vi <- 1/data_cc$y1 + 1/data_cc$c1

lnOR_nn <- NN_LMM(
  yi=yi, vi=vi, 
  mu.upper = 5, tau.upper = 3,
  init.vals = c(-3,0.1))
lnOR_nn_lb <- lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_nn$mu[2]
lnOR_nn_ub <- lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_nn$mu[2]

sprintf("theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_nn$mu[1], lnOR_nn_lb, lnOR_nn_ub, lnOR_nn$mu[2])
sprintf("tau (SE): %.3f (%.3f)", 
        lnOR_nn$tau[1], lnOR_nn$tau[2])

## SENSITIVITY ANALYSIS FOR PUBLICATION BIAS IN META-ANALYSIS-------------------
## (Pnmax = Psemin = p)---------------------------------------------------------
## 
## PROPOSAL BASED ON BN-GLMM
p_sa <- seq(0.9, 0.1, -0.1)
lnOR_copas_BNGLMM <- vapply(
  p_sa, 
  function(p) {
    mod <- copas_BNGLMM_prop(y1=y1, n1=n1, Pnmax = 0.99, Pnmin = p, 
                        rho.init = 0.1, mu.upper = 6, tau.upper = 3, integ.limit = 20,
                        init.vals = c(lnOR_bn$mu[1], lnOR_bn$tau[1]))
    mu_lb <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$a, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "a0"=0, "a1"=0, "converge"=0))

colnames(lnOR_copas_BNGLMM) <- paste0("p = ", p_sa)
lnOR_copas_BNGLMM

## COPAS-HECKMAN-TYPE SELECTION FUNCTION
lnOR_copas_NNLMM <- vapply(
  p_sa, 
  function(p) {
    mod <- copas_NNLMM(yi = yi, vi = vi, Psemax = p, Psemin = 0.99, 
                       rho.init = 0.1, mu.upper = 6, tau.upper = 3,
                       init.vals = c(lnOR_nn$mu[1], lnOR_nn$tau[1]))
    mu_lb <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))

colnames(lnOR_copas_NNLMM) <- paste0("p = ", p_sa)
lnOR_copas_NNLMM

## NUMBER OF THE UNPUBLISHED
M_propos <- sapply(p_sa, function(p) {
  
  P_max <- p
  P_min <- 0.99

  n_min <- min(n1) 
  n_max <- max(n1)
  
  a1 <- (qnorm(P_max)-qnorm(P_min))/(sqrt(n_max)-sqrt(n_min))
  a0 <- qnorm(P_max)-a1*sqrt(n_max)
  sum((1 - pnorm(a0+a1*sqrt(n1)))/pnorm(a0+a1*sqrt(n1))) 
  
})


M_copas <- sapply(p_sa, function(p) {
  
  P_max <- 0.99
  P_min <- p

  se_min_inv <- 1/sqrt(min(vi)) 
  se_max_inv <- 1/sqrt(max(vi))
  
  gamma1 <- (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 <- qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi)))/pnorm(gamma0+gamma1/sqrt(vi)))
  
})


## Table 1
tab1 <- data.frame(
  BN = t(lnOR_copas_BNGLMM[c(1,3,4),]),
  NN = t(lnOR_copas_NNLMM[c(1,3,4),]),
  M.c = round(M_copas), M.p = round(M_propos)
)

tab1_p1 <- c(
  BN = lnOR_bn$mu[1], BN.mu.lb = unname(lnOR_bn_lb), BN.mu.ub = unname(lnOR_bn_ub),
  NN = lnOR_nn$mu[1], NN.mu.lb = unname(lnOR_nn_lb), NN.mu.ub = unname(lnOR_nn_ub),
  M.c = 0, M.p = 0
)

tab1_all <- rbind("p = 1" = tab1_p1, tab1)
tab1_all$pnmax <- c(NA, p_sa)
tab1_all$pnmin <- c(NA, rep(0.99, 9))

## SAVE RESULTS1
## 
# save(lnOR_copas_BNGLMM, lnOR_copas_NNLMM,
#      M_propos, M_copas, tab1_all,
#      lnOR_bn, lnOR_nn,
#      file = "example-supp-bias1.RData")


## (Pnmin = Psemax = p)---------------------------------------------------------
## 
## PROPOSAL BASED ON BN-GLMM
lnOR_copas_BNGLMM <- vapply(
  p_sa, 
  function(p) {
    mod <- copas_BNGLMM_prop(y1=y1, n1=n1, Pnmax = 0.99, Pnmin = p, 
                        rho.init = 0.1, mu.upper = 6, tau.upper = 3, integ.limit = 20,
                        init.vals = c(lnOR_bn$mu[1], lnOR_bn$tau[1]))
    mu_lb <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$a, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "a0"=0, "a1"=0, "converge"=0))

colnames(lnOR_copas_BNGLMM) <- paste0("p = ", p_sa)
lnOR_copas_BNGLMM

## COPAS-HECKMAN-TYPE SELECTION FUNCTION
lnOR_copas_NNLMM <- vapply(
  p_sa, 
  function(p) {
    mod <- copas_NNLMM(yi = yi, vi = vi, Psemax = p, Psemin = 0.99, 
                       rho.init = 0.1, mu.upper = 6, tau.upper = 3,
                       init.vals = c(lnOR_nn$mu[1], lnOR_nn$tau[1]))
    mu_lb <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub <- mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))

colnames(lnOR_copas_NNLMM) <- paste0("p = ", p_sa)
lnOR_copas_NNLMM

## NUMBER OF THE UNPUBLISHED
M_propos <- sapply(p_sa, function(p) {
  
  P_max <- 0.99
  P_min <- p
  
  n_min <- min(n1) 
  n_max <- max(n1)
  
  a1 <- (qnorm(P_max)-qnorm(P_min))/(sqrt(n_max)-sqrt(n_min))
  a0 <- qnorm(P_max)-a1*sqrt(n_max)
  sum((1 - pnorm(a0+a1*sqrt(n1)))/pnorm(a0+a1*sqrt(n1))) 
  
})


M_copas <- sapply(p_sa, function(p) {
  
  P_max <- p
  P_min <- 0.99

  se_min_inv <- 1/sqrt(min(vi)) 
  se_max_inv <- 1/sqrt(max(vi))
  
  gamma1 <- (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 <- qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi)))/pnorm(gamma0+gamma1/sqrt(vi)))
  
})


## Table 2
tab2 <- data.frame(
  BN = t(lnOR_copas_BNGLMM[c(1,3,4),]),
  NN = t(lnOR_copas_NNLMM[c(1,3,4),]),
  M.c = round(M_copas), M.p = round(M_propos)
)


tab2_all <- rbind("p = 1" = tab1_p1, tab2)
tab2_all$pnmin <- c(NA, p_sa)
tab2_all$pnmax <- c(NA, rep(0.99, 9))

## SAVE RESULTS2
## 
save(lnOR_copas_BNGLMM, lnOR_copas_NNLMM,
     M_propos, M_copas, tab2_all,
     lnOR_bn, lnOR_nn,
     file = "example-supp-bias2.RData")

