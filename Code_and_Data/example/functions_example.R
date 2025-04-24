library(extraDistr)
library(nleqslv)

#' General purpose error handling function
#'
#' @param e error output of the tryCatch function
#' @return error message in the console; list of missing values
errorfun <- function(e) {
  message("Error occurs in the code. Here's the original error message:")
  message(e)
  # Choose a return value in case of error
  return(list(b = NA, se = NA, ci.lb = NA, ci.ub = NA, tau2 = NA))
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
