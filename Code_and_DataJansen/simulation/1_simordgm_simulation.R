#---------------------------------------------------------------------------------------------
# SIMULATION
#---------------------------------------------------------------------------------------------

# this script runs the simulation for all conditions and saves the simulation results of each
# condition in a file named "simor_i.rds", where i is the index of the simulation condition
# the "simor_i.rds" files are saved in a folder named "results", which will be created
# if it does not exist already.
# processing of these files is carried out by the script "2_simordgm_processing.R"

# load required packages -----------------------------------------------------------------
library(doParallel) # package for parallel computing
library(doRNG) # makes parallel simulation reproducible
library(tidyr) # data reformatting
library(dplyr) # data reformatting
library(lme4) # glmms
library(metafor) # iv model, nchg model
library(BiasedUrn) # required to fit noncentral hypergeometric-normal model using metafor
library(extraDistr) # for beta distribution
library(nleqslv) # non-linear equation solver

# set working directory to folder which contains the file "functions_simordgm.R"
setwd("~/simulation")

source("functions_simordgm.R")

# create results directory if it does not exist already
dir.create("results")

# set up parallel computing on cluster -------------------------------------------------------
# note that this must be adapted if you run the code on your own computer
# for instance, replace "as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE"))" by "detectCores()-1"
cores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE"))
cluster <- makeCluster(cores)
registerDoParallel(cluster, cores = cores)

# save session Info -------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo_SIMORdgm.txt")

# set simulation parameters ------------------------------------------------------------------
design <- expand.grid(
  true_or = c(0.5, 1, 2),
  p_control = c(0.01, 0.05),
  n_studies = c(5, 10, 30),
  n_sample = c(50, 100, 200, 500),
  group_ratio = c(1, 2),
  tau = c(0, 0.3, 0.6, 1),
  dgm = c("pCFixed", "pRandom")
) %>% arrange(p_control, n_studies, n_sample)

save(design, file = "./results/design_simordgm.RData")

n_conditions <- nrow(design)
n_simulations <- 1000

# simulation --------------------------------------------------------------------------------

I <- 1:n_conditions
for (i in I) {

  # print where we are
  message("We are in the ", i, "th condition of ", n_conditions, ".")

  registerDoRNG(145) # set seed to make simulation reproducible

  # loop through simulation trials within each condition
  # this is parallelized
  # results for the ith condition for all trials will be saved in results folder
  results_simulation_trials <- foreach(j = 1:n_simulations, .errorhandling = "pass") %dopar% {
    require(dplyr)
    require(tidyr)
    require(lme4)
    require(metafor)
    require(BiasedUrn)
    require(extraDistr)
    require(nleqslv)
    require(numDeriv)

    # simulate data set for the trial --------------------------------------------------------

    data <- generate_data(
      theta = log(design$true_or[i]),
      tau = design$tau[i],
      p0 = design$p_control[i],
      n = design$n_sample[i],
      gr = design$group_ratio[i],
      k = design$n_studies[i],
      dgm = design$dgm[i]
    )


    # reformat data ------------------------------------------------------------------------------

    # data iv model

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


    # estimate models ----------------------------------------------------------------------------

    # initiate model list
    model_list <- list()

    # IV models ------------------------

    start <- Sys.time()
    m_sj <- tryCatch(
      {
        rma(yi = yi, vi = vi, data = data_cc, method = "SJ")
      },
      error = errorfun
    )
    time_sj <- Sys.time() - start

    model_list[["m_sj"]] <- list(model = m_sj, runtime = time_sj)


    # IV model with KNHA and SJ estimator ------------------------------------------

    start <- Sys.time()
    m_sj_knha <- tryCatch(
      {
        rma(yi = yi, vi = vi, data = data_cc, method = "SJ", test = "knha")
      },
      error = errorfun
    )
    time_sj_knha <- Sys.time() - start

    model_list[["m_sj_knha"]] <- list(model = m_sj_knha, runtime = time_sj_knha)


    # Mantel-Haenszel model --------------------------------------------------------

    start <- Sys.time()
    m_mh <- tryCatch(
      {
        rma.mh(
          ai = event_treat, bi = noevent_treat,
          ci = event_control, di = noevent_control,
          data = data, to = "none"
        )
      },
      error = errorfun
    )
    time_mh <- Sys.time() - start

    model_list[["m_mh"]] <- list(model = m_mh, runtime = time_mh)



    # Binomial models ------------------

    # binomial-normal model (fixed intercept, usual parametrization)
    start <- Sys.time()
    m_binom_fi <- tryCatch(
      {
        glmer(cbind(event_count, n - event_count) ~ factor(study) + factor(dummy_group) +
          (0 + dummy_group | study),
        data = data_long,
        family = binomial(link = "logit")
        )
      },
      error = errorfun
    )
    time_binom_fi <- Sys.time() - start

    # save model
    model_list[["m_binom_fi"]] <- list(model = m_binom_fi, runtime = time_binom_fi)


    # binomial-normal model (fixed intercept, alternative parametrization)
    start <- Sys.time()
    m_binom_fi2 <- tryCatch(
      {
        glmer(cbind(event_count, n - event_count) ~ factor(study) + factor(dummy_group) +
          (0 + dummy_group_alt | study),
        data = data_long,
        family = binomial(link = "logit")
        )
      },
      error = errorfun
    )
    time_binom_fi2 <- Sys.time() - start

    # save models
    model_list[["m_binom_fi2"]] <- list(model = m_binom_fi2, runtime = time_binom_fi2)


    # random intercept models ----------------

    # binomial-normal model (random intercept, usual parametrization)
    start <- Sys.time()
    m_binom_ri <- tryCatch(
      {
        glmer(cbind(event_count, n - event_count) ~ factor(dummy_group) + (dummy_group || study),
          data = data_long,
          family = binomial(link = "logit")
        )
      },
      error = errorfun
    )
    time_binom_ri <- Sys.time() - start

    # save model
    model_list[["m_binom_ri"]] <- list(model = m_binom_ri, runtime = time_binom_ri)

    # binomial-normal model (random intercept, alternative parametrization)
    start <- Sys.time()
    m_binom_ri2 <- tryCatch(
      {
        glmer(cbind(event_count, n - event_count) ~ factor(dummy_group) + (dummy_group_alt || study),
          data = data_long, family = binomial(link = "logit")
        )
      },
      error = errorfun
    )
    time_binom_ri2 <- Sys.time() - start

    # save model
    model_list[["m_binom_ri2"]] <- list(model = m_binom_ri2, runtime = time_binom_ri2)


    # Beta-binomial model ------------
    theta_init_beta <- c(-1, 0, 1) # starting values
    names(theta_init_beta) <- c("b0", "b1", "prec")

    start <- Sys.time()
    m_betabinom <- tryCatch(
      {
        optim(
          par = theta_init_beta,
          fn = likelihood_betabin,
          event_treat = data$event_treat,
          noevent_treat = data$noevent_treat,
          event_control = data$event_control,
          noevent_control = data$noevent_control,
          hessian = TRUE
        )
      },
      error = errorfun
    )
    time_betabinom <- Sys.time() - start

    model_list[["m_betabinom"]] <- list(model = m_betabinom, runtime = time_betabinom)


    # Nchg model ---------------------

    start <- Sys.time()
    m_nchg <- tryCatch(
      {
        rma.glmm(
          ai = event_treat, bi = noevent_treat,
          ci = event_control, di = noevent_control,
          measure = "OR",
          data = data,
          model = "CM.EL",
          control = list(
            optCtrl = list(
              maxit = 20000, reltol = 0.0001
            ), optmethod = "Nelder-Mead",
            dnchgcalc = "dFNCHypergeo"
          )
        )
      },
      error = errorfun
    )
    time_nchg <- Sys.time() - start

    model_list[["m_nchg"]] <- list(model = m_nchg, runtime = time_nchg)

    # Conditional binomial model ---------------------------------------------------
    start <- Sys.time()
    m_condbinom <- tryCatch(
      {
        rma.glmm(
          measure = "OR", ai = event_treat, bi = noevent_treat,
          ci = event_control, di = noevent_control,
          data = data, model = "CM.AL", method = "ML"
        )
      },
      error = errorfun
    )
    time_condbinom <- Sys.time() - start

    model_list[["m_condbinom"]] <- list(model = m_condbinom, runtime = time_condbinom)

    # return model list and data for j-th simulation replication
    return(list(model_list = model_list, data = data))
  } # end loop simulation trial

  # save results of simulation trial in results folder
  save_to_file <- paste0("./results/simor_", i, ".rds")
  saveRDS(results_simulation_trials, file = save_to_file)
} # end loop condition


stopCluster(cluster)
