#---------------------------------------------------------------------------------------------
# 5. EXAMPLES: Thomas et al. 2015
#---------------------------------------------------------------------------------------------

library("lme4")
library("metafor")
library("tidyverse")
library("xtable")

setwd("~/example")

source("functions_example.R")

# Prepare data -------------------------------------------------------------------------------
load("data_thomas.RData")

# 1 Data Thomas ------------------------------------------------------------------------------
outcomes <- c("Suicidal Ideation", "Aggression")
n_outcomes <- length(outcomes)

for (i in 1:n_outcomes) {
  tryCatch(
    {
      data <- subset(data_thomas, outcome == outcomes[i])

      # standard continuity correction (+ 0.5 added to single-zero studies)
      data_cc <- escalc("OR",
        ai = event_treat, bi = noevent_treat,
        ci = event_control, di = noevent_control,
        data = data, add = 0.5, to = "only0", drop00 = TRUE
      ) %>%
        filter(!is.na(yi)) # deletes double 0 studies

      # data long format
      data_long <- data %>%
        gather(c(event_treat, event_control), key = "group", value = "event_count") %>%
        mutate(
          group = as.factor(group),
          group = recode_factor(group, event_control = "control", event_treat = "treat"),
          dummy_group = ifelse(group == "treat", 1, 0),
          dummy_group_rev = ifelse(group == "treat", 0, 1),
          dummy_group_alt = ifelse(group == "treat", 0.5, -0.5),
          n = ifelse(group == "treat", n_treat, n_control)
        ) %>%
        select(study, group, event_count, dummy_group, dummy_group_rev, dummy_group_alt, n)

      # IV RE model ----------------------------------------------------------------------------------

      # SJ:
      m_sj <- rma(yi = yi, vi = vi, data = data_cc, method = "SJ")

      # SJ-HKSJ:
      m_sj_hksj <- rma(yi = yi, vi = vi, data = data_cc, method = "SJ", test = "knha")

      # BinFI models ---------------------------------------------------------------------------------

      # BinFI:
      m_binom_fi <- glmer(cbind(event_count, n - event_count) ~ factor(study) + factor(dummy_group) +
        (0 + dummy_group | study),
      data = data_long,
      family = binomial(link = "logit")
      )

      CI_binom_fi <- confint(m_binom_fi, method = "Wald")["factor(dummy_group)1", ]


      # BinFI2:

      m_binom_fi2 <- glmer(cbind(event_count, n - event_count) ~ factor(study) + factor(dummy_group) +
        (0 + dummy_group_alt | study),
      data = data_long,
      family = binomial(link = "logit")
      )

      CI_binom_fi2 <- confint(m_binom_fi2, method = "Wald")["factor(dummy_group)1", ]


      # BinRI models ---------------------------------------------------------------------------------

      # BinRI:
      m_binom_ri <- glmer(cbind(event_count, n - event_count) ~ factor(dummy_group) + (dummy_group || study),
        data = data_long,
        family = binomial(link = "logit")
      )

      CI_binom_ri <- confint(m_binom_ri, method = "Wald")["factor(dummy_group)1", ]

      # BinRI2:
      m_binom_ri2 <- glmer(cbind(event_count, n - event_count) ~ factor(dummy_group) + (dummy_group_alt || study),
        data = data_long, family = binomial(link = "logit")
      )

      CI_binom_ri2 <- confint(m_binom_ri2, method = "Wald")["factor(dummy_group)1", ]


      # Betabin model --------------------------------------------------------------------------------

      theta_init_beta <- c(-1, 0, 1)
      names(theta_init_beta) <- c("b0", "b1", "prec")

      m_betabin <- optim(
        par = theta_init_beta,
        fn = likelihood_betabin,
        event_treat = data$event_treat,
        noevent_treat = data$noevent_treat,
        event_control = data$event_control,
        noevent_control = data$noevent_control,
        hessian = TRUE
      )

      se_betabin <- sqrt(solve(m_betabin$hessian)["b1", "b1"])

      CI_betabin <- c(m_betabin$par["b1"] - qnorm(0.975) * se_betabin, m_betabin$par["b1"] + qnorm(0.975) * se_betabin)

      # Nchg model -----------------------------------------------------------------------------------

      m_nchg <- rma.glmm(
        ai = event_treat, bi = noevent_treat,
        ci = event_control, di = noevent_control,
        measure = "OR",
        data = data,
        model = "CM.EL",
        control = list(optCtrl = list(
          maxit = 20000, reltol = 0.0001
        ), 
          optmethod = "Nelder-Mead", 
          dnchgcalc = "dFNCHypergeo")
      )

      # Condbin model ---------------------------------------------------------------------------------

      m_condbinom <- rma.glmm(
        measure = "OR", ai = event_treat, bi = noevent_treat,
        ci = event_control, di = noevent_control,
        data = data, model = "CM.AL", method = "ML",
        control = list(optimizer = "optim")
      )

      # Mantel-Haenszel model -------------------------------------------------------------------------

      m_mh <- rma.mh(
        ai = event_treat, bi = noevent_treat,
        ci = event_control, di = noevent_control,
        data = data, to = "none"
      )

      # Original Analysis: Peto OR --------------------------------------------------------------------
      m_peto <- rma.peto(ai = data$event_treat, bi = data$noevent_treat, ci = data$event_control, di = data$noevent_control)


      # Summarize Results in Table using xtable -------------------------------------------------------

      models <- c(
        "Peto OR", "IV_SJ", "IV_SJ_HKSJ",
        "BinFI", "BinFI2", "BinRI", "BinRI2",
        "Nchg", "Condbin",
        "Betabin", "MH"
      )

      logORs <- c(
        m_peto$b, m_sj$b, m_sj_hksj$b,
        summary(m_binom_fi)$coefficients["factor(dummy_group)1", 1],
        summary(m_binom_fi2)$coefficients["factor(dummy_group)1", 1],
        summary(m_binom_ri)$coefficients["factor(dummy_group)1", 1],
        summary(m_binom_ri2)$coefficients["factor(dummy_group)1", 1],
        m_nchg$b, m_condbinom$b,
        as.numeric(m_betabin$par["b1"]), m_mh$b
      )
      logORs <- round(logORs, 2)

      CI_lower <- c(
        m_peto$ci.lb, m_sj$ci.lb, m_sj_hksj$ci.lb,
        as.numeric(CI_binom_fi[1]),
        as.numeric(CI_binom_fi2[1]),
        as.numeric(CI_binom_ri[1]),
        as.numeric(CI_binom_ri2[1]),
        m_nchg$ci.lb, m_condbinom$ci.lb,
        as.numeric(CI_betabin[1]), m_mh$ci.lb
      )
      CI_lower <- round(CI_lower, 2)

      CI_upper <- c(
        m_peto$ci.ub, m_sj$ci.ub, m_sj_hksj$ci.ub,
        as.numeric(CI_binom_fi[2]),
        as.numeric(CI_binom_fi2[2]),
        as.numeric(CI_binom_ri[2]),
        as.numeric(CI_binom_ri2[2]),
        m_nchg$ci.ub, m_condbinom$ci.ub,
        as.numeric(CI_betabin[2]), m_mh$ci.ub
      )
      CI_upper <- round(CI_upper, 2)

      CIs <- paste0("[", CI_lower, ", ", CI_upper, "]", sep = "")

      tau2s <- c(
        NA, m_sj$tau2, m_sj_hksj$tau2, 
        VarCorr(m_binom_fi)$study["dummy_group", "dummy_group"]^2,
        VarCorr(m_binom_fi2)$study["dummy_group_alt", "dummy_group_alt"]^2,
        VarCorr(m_binom_ri)$study.1["dummy_group", "dummy_group"]^2,
        VarCorr(m_binom_ri2)$study.1["dummy_group_alt", "dummy_group_alt"]^2,
        m_nchg$tau2, m_condbinom$tau2,
        NA, NA
      )
      tau2s <- round(tau2s, 2)


      example_df <- data.frame(
        Data = rep("Thomas et al. 2015", 11),
        Outcome = rep(outcomes[i], 11),
        Model = models,
        logOR = logORs,
        CI = CIs,
        tau2 = tau2s
      )

      caption <- paste0("Results for the Data from Thomas et al. 2015 for the Outcome ", outcomes[i], sep = "")
      texfile_name <- paste0("tab_thomas2015_", outcomes[i], ".tex", sep = "")
      names(example_df)[c(4, 6)] <- c("$log(\\widehat{OR})$", "$\\hat{\\tau}^2$")

      example_xtable <- xtable(example_df[3:6], caption = caption)
      print(xtable(example_df[3:6], type = "latex", caption = caption), sanitize.text.function = function(x) {
        x
      }, caption.placement = "top", include.rownames = FALSE, file = texfile_name)


      # file_name <- paste0("thomas2015_reanalysis_", i, ".RData")
      # save(example_df, file = file_name)
    },
    error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
}
