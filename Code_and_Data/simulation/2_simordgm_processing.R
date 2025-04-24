#---------------------------------------------------------------------------------------------
# DATA PROCESSING
#---------------------------------------------------------------------------------------------

# this script processes the simulation results.
# it loads the "simor_i.rds" files from the results folder
# and produces a data frame which includes the results of interest
# this data frame is saved in the results folder as "simor_summary.RData".
# it is required to (re)produce the figures from the manuscript and supplement, which can be
# done using the scripts "3_simordgm_plots_manuscript.R" and "4_simordgm_plots_supplement.R"

library(lme4)
library(dplyr)
library(tidyr)

# set working directory to folder which contains simor_res_i.rds files AND design_simordgm.Rdata file
setwd("~/simulation/results")

writeLines(capture.output(sessionInfo()), "sessionInfo_SIMORdgm_processing.txt")


# import design file
load("design_simordgm.RData")

# get number of conditions
n_cond_design <- nrow(design)

# list result files
res_files <- list.files(path = ".", pattern = "simor_")

# get numbers of simulated conditions
sim_conditions <- sort(as.integer(sub(".*_(\\d+).*", "\\1", res_files)))

n_cond <- length(sim_conditions)

# number of missing conditions
missing_cond <- n_cond_design - n_cond

# warn if condition files are missing:
if (missing_cond > 0) {
  warning(paste(missing_cond, "condition files are missing"))
}

# model params
params <- c("error", "logOR", "se.logOR", "ci.lb", "ci.ub", "runtime")

# DATA EVALUATION -----------------------------------------------------------------------------------

# loop through conditions:
results_full <- lapply(
  sim_conditions,
  function(i) {
    message("Evaluating condition ", i)

    file_name <- paste0("simor_", i, ".rds")


    # load result file for i-th condition
    res_i <- readRDS(file_name)

    n_sim_i <- length(res_i) # no of trials which were simulated
    models_i <- c("m_sj", "m_sj_knha", "m_mh", "m_nchg", "m_condbinom",
                  "m_binom_fi", "m_binom_fi2", "m_binom_ri", "m_binom_ri2",
                  "m_betabinom"
                  )

    # loop through simulation trials:
    results_j <- lapply(
      1:n_sim_i,
      function(j) {
        message("Evaluating trial ", j)

        # data evaluation
        data <- res_i[[j]]$data

        n_dz <- sum(data$event_treat == 0 & data$event_control == 0)
        all_dz <- (n_dz == nrow(data))

        n_sz_treat <- sum(data$event_treat == 0 & !(data$event_control == 0))
        all_sz_treat <- (sum(data$event_treat == 0) == nrow(data))

        n_sz_control <- sum(data$event_control == 0 & !(data$event_treat == 0))
        all_sz_control <- (sum(data$event_control == 0) == nrow(data))

        total_treat <- sum(data$event_treat)
        total_control <- sum(data$event_control)


        results_j <- data.frame(
          condition = i,
          sim_trial = j,
          true_logOR = log(design$true_or[i]),
          true_tau = design$tau[i],
          p_control = design$p_control[i],
          n_studies = design$n_studies[i],
          n_sample = design$n_sample[i],
          group_ratio = design$group_ratio[i],
          dgm = design$dgm[i],
          n_dz = n_dz, all_dz = all_dz,
          n_sz_treat = n_sz_treat, all_sz_treat = all_sz_treat,
          n_sz_control = n_sz_control, all_sz_control = all_sz_control,
          total_treat = total_treat, total_control = total_control
        )

        # model evaluation

        # # loop through models:
        model_params_j <- for (m in 1:length(models_i)) {
          model <- models_i[m]
          model_res <- res_i[[j]]$model_list[[model]]$model

          param_names <- paste(model, params, sep = "_")
          results_j[, param_names] <- NA

          results_j[[paste(model, "runtime", sep = "_")]] <- as.numeric(res_i[[j]]$model_list[[model]]$runtime, units = "secs")

          # IV models, nchg model
          if (model %in% c("m_sj", "m_sj_knha", "m_mh", "m_nchg", "m_condbinom")) {
            model_error <- ifelse(!is.null(model_res[["error"]]), as.character(model_res$error), "None")
            results_j[[paste(model, "error", sep = "_")]] <- model_error

            if (model_error == "None") {
              results_j[[paste(model, "logOR", sep = "_")]] <- as.numeric(model_res$beta)
              results_j[[paste(model, "se.logOR", sep = "_")]] <- model_res$se
              results_j[[paste(model, "ci.lb", sep = "_")]] <- model_res$ci.lb
              results_j[[paste(model, "ci.ub", sep = "_")]] <- model_res$ci.ub
            } # model error end
          } # rma model end

          # Binomial models
          if (model %in% c("m_binom_fi", "m_binom_fi2", "m_binom_ri", "m_binom_ri2")) {
            model_error <- ifelse(class(model_res) == "list", as.character(model_res$error), "None")
            results_j[[paste(model, "error", sep = "_")]] <- model_error

            if (model_error == "None") {
              results_j[[paste(model, "logOR", sep = "_")]] <- as.numeric(fixef(model_res)["factor(dummy_group)1"])
              results_j[[paste(model, "se.logOR", sep = "_")]] <- tryCatch(
                {
                  coef(summary(model_res))["factor(dummy_group)1", "Std. Error"]
                },
                error = function(cond) {
                  return(NA)
                }
              )
              results_j[[paste(model, "ci.lb", sep = "_")]] <- results_j[[paste(model, "logOR", sep = "_")]] - qnorm(0.975) * results_j[[paste(model, "se.logOR", sep = "_")]]
              results_j[[paste(model, "ci.ub", sep = "_")]] <- results_j[[paste(model, "logOR", sep = "_")]] + qnorm(0.975) * results_j[[paste(model, "se.logOR", sep = "_")]]
            } # model error end
          } # binomial model end

          # Betabinomial model
          if (model %in% c("m_betabinom")) {
            model_error <- ifelse(!is.null(model_res[["error"]]), as.character(model_res$error), "None")
            results_j[[paste(model, "error", sep = "_")]] <- model_error

            if (model_error == "None") {
              results_j[[paste(model, "logOR", sep = "_")]] <- as.numeric(model_res$par["b1"])
              results_j[[paste(model, "se.logOR", sep = "_")]] <- tryCatch(
                {
                  sqrt(solve(model_res$hessian)["b1", "b1"])
                },
                error = function(cond) {
                  return(NA)
                }
              )
              results_j[[paste(model, "ci.lb", sep = "_")]] <- as.numeric(model_res$par["b1"]) - qnorm(0.975) * results_j[[paste(model, "se.logOR", sep = "_")]]
              results_j[[paste(model, "ci.ub", sep = "_")]] <- as.numeric(model_res$par["b1"]) + qnorm(0.975) * results_j[[paste(model, "se.logOR", sep = "_")]]
            } # model_error end
          } # betabinomial model end

        } # model loop end



        return(results_j)
      } # lapply inner function end
    ) # lapply inner end

    results_j_df <- do.call(rbind, results_j)

    # save intermediate results for condition i:
    results_df_name <- paste0("results_", i, ".RData")
    save(results_j_df, file = results_df_name)

    # create summary which is used to make figures for manuscript and supplement:
    results_j_summary <- results_j_df %>%
      mutate(
        sum_all_dz = sum(all_dz),
        sum_all_sz_treat = sum(all_sz_treat),
        sum_all_sz_control = sum(all_sz_control),
        mean_n_dz = mean(n_dz),
        mean_n_sz_treat = mean(n_sz_treat),
        mean_n_sz_control = mean(n_sz_control),
        mean_total_treat = mean(total_treat),
        mean_total_control = mean(total_control)
      ) %>%
      pivot_longer(
        cols = m_sj_error:m_betabinom_runtime, names_to = c("model", ".value"),
        names_pattern = "m_(.*)_(.*)"
      ) %>%
      group_by(model) %>%
      mutate(
        estimability_not_total0 = sum(!is.na(se.logOR) & all_dz == FALSE) / sum(all_dz == FALSE),
        estimability_not_total0sz = sum(!is.na(se.logOR) & all_dz == FALSE & all_sz_treat == FALSE & all_sz_control == FALSE) / sum(all_dz == FALSE & all_sz_treat == FALSE & all_sz_control == FALSE),
        error = mean(error != "None"),
        estimability = mean(!is.na(se.logOR))
      ) %>%
      filter(!is.na(se.logOR)) %>%
      mutate(
        ci.lb = ifelse(is.na(ci.lb), logOR - qnorm(0.975) * se.logOR, ci.lb),
        ci.ub = ifelse(is.na(ci.ub), logOR + qnorm(0.975) * se.logOR, ci.ub)
      ) %>%
      group_by(model) %>%
      summarise(
        condition = i,
        true_logOR = true_logOR[1],
        p_control = p_control[1],
        true_tau = true_tau[1],
        n_studies = n_studies[1],
        n_sample = n_sample[1],
        group_ratio = group_ratio[1],
        dgm = dgm[1],
        sum_all_dz = sum_all_dz[1],
        sum_all_sz_treat = sum_all_sz_treat[1],
        sum_all_sz_control = sum_all_sz_control[1],
        mean_n_dz = mean_n_dz[1],
        mean_n_sz_treat = mean_n_sz_treat[1],
        mean_n_sz_control = mean_n_sz_control[1],
        mean_total_treat = mean_total_treat[1],
        mean_total_control = mean_total_control[1],
        error = error[1],
        estimability = estimability[1],
        estimability_not_total0 = estimability_not_total0[1],
        estimability_not_total0sz = estimability_not_total0sz[1],
        mean_bias_logOR = mean(logOR - true_logOR[1], na.rm = TRUE),
        median_bias_logOR = median(logOR - true_logOR[1], na.rm = TRUE),
        var_logOR = var(logOR, na.rm = TRUE),
        mcse_logOR = sqrt(var(logOR, na.rm = TRUE)/(estimability[1]*1000)),
        mse_logOR = mean((logOR - true_logOR[1])^2, na.rm = TRUE),
        coverage_logOR = mean((true_logOR[1] <= ci.ub) & (true_logOR[1] >= ci.lb), na.rm = TRUE),
        mean_ciwidth_logOR = mean(ci.ub - ci.lb, na.rm = TRUE),
        median_ciwidth_logOR = median(ci.ub - ci.lb, na.rm = TRUE),
        runtime = mean(runtime)
      )

    return(results_j_summary)
  } # lapply outer function end
) # lapply outer end

results_full <- do.call(rbind, results_full)
save(results_full, file = "simor_summary.RData")
