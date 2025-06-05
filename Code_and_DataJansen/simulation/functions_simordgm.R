library(extraDistr)
library(nleqslv)

#' General purpose error handling function
#'
#' @param cond error output of the tryCatch function
#' @return returns a list of the errors
errorfun <- function(cond) {
  # Choose a return value in case of error
  return(list(error = cond))
}


#' Generate data using the algorithm pCFixed or pRandom
#'
#' @param theta true log OR
#' @param tau true between-study standard deviation
#' @param p0 true average control group probability
#' @param n total sample size (ntreat + ncontrol)
#' @param gr group ratio (ntreat/ncontrol)
#' @param k number of studies
#' @param dgm data-generating model; options: "pCFixed" or "pRandom"
#'
#' @return data frame containing the simulated data
generate_data <- function(theta, tau, p0 = NA, n, gr, k, dgm = "pRandom") {

  # generate sample sizes
  n0i <- rbinom(k, n, 1 / (1 + gr)) # decides for each observation
  n1i <- n - n0i


  if (dgm == "pRandom") {
    p1_init <- p0 * exp(theta) / ( 1 - p0 + p0 * exp(theta))

    mu0 <- log(p0 / (1 - p0))
    mu1 <- log(p1_init / (1 - p1_init))

    logit0i <- rnorm(k, mu0, tau / sqrt(2))
    logit1i <- rnorm(k, mu1, tau / sqrt(2))

    p0i <- 1 / (1 + exp(-logit0i))
    p1i <- 1 / (1 + exp(-logit1i))

    y0i <- rbinom(k, n0i, p0i)
    y1i <- rbinom(k, n1i, p1i)
  }


  if (dgm == "pCFixed") {
    thetai <- rnorm(k, mean = theta, sd = tau)

    p0i <- p0 # (held constant at the moment)
    p1i <- p0i * exp(thetai) / (1 - p0i + p0i * exp(thetai))

    y0i <- rbinom(k, n0i, p0i)
    y1i <- rbinom(k, n1i, p1i)
  }


  data <- data.frame(
    study = 1:k,
    event_treat = y1i,
    noevent_treat = n1i - y1i,
    event_control = y0i,
    noevent_control = n0i - y0i,
    n_treat = n1i,
    n_control = n0i
  )

  attr(data, "theta") <- theta
  attr(data, "tau") <- tau
  attr(data, "k") <- k
  attr(data, "n") <- n
  attr(data, "gr") <- gr
  attr(data, "dgm") <- dgm

  data
}


#' Likelihood of the beta-binomial model
#'
#' @param theta named vector of starting values (with names "b0", "b1" and "prec")
#' @param event_treat vector which contains the numbers of events in the treatment
#' group for all k studies
#' @param noevent_treat vector which contains the numbers of non-events in the
#' treatment group for all k studies
#' @param event_control vector which contains the numbers of events in the control
#' group for all k studies
#' @param noevent_control vector which contains the numbers of non-events in the
#' control group for all k studies
#' @return returns the value of the beta-binomial likelihood for parameters theta
#' and data specified by event_treat, noevent_treat, event_control and noevent_control
likelihood_betabin <- function(theta,
                                 event_treat, noevent_treat,
                                 event_control, noevent_control) {
  b0 <- theta[1]
  b1 <- theta[2]
  prec <- theta[3]
  total_treat <- event_treat + noevent_treat
  total_control <- event_control + noevent_control
  mu_treat <- exp(b0 + b1) / (1 + exp(b0 + b1))
  mu_control <- exp(b0) / (1 + exp(b0))

  # mean precision parameterization of beta-binomial distribution
  alpha_treat <- mu_treat * prec
  beta_treat <- (1 - mu_treat) * prec
  alpha_control <- mu_control * prec
  beta_control <- (1 - mu_control) * prec

  to_sum_treat <- dbbinom(
    event_treat,
    size = total_treat, alpha = alpha_treat,
    beta = beta_treat, log = TRUE
  )
  to_sum_control <- dbbinom(
    event_control,
    size = total_control, alpha = alpha_control,
    beta = beta_control, log = TRUE
  )

  sum_total <- sum(to_sum_treat) + sum(to_sum_control)

  return(-1 * sum_total)
}
