#---------------------------------------------------------------------------------------------
# PLOTS (SUPPLEMENT)
#---------------------------------------------------------------------------------------------

# this script produces the plots for the supplement.
# it requires that the file "simordgm_summary.RData", which contains the summarized
# simulation results, exists in the results folder in the working directory.
# this script creates a directory called plots_supp in which the
# figures are saved as pdf files.

library(tidyverse)
library(ggh4x) # nested plots
library(ggthemes)
library(scales) # required to manually adjust the scales of the barplots

# set working directory to folder which contains the folder "results"
# the results folder has to contain at least the files "results_simordgm.RData" and "design_simordgm.RData"
setwd("~/results")

# load data and design file
load("simor_summary.RData")

# create  directory for plots if it does not exit already
dir.create("plots_supp")
plot_dir <- "./plots_supp"

# round value of the true log OR to two decimals
results_full$true_logOR <- round(results_full$true_logOR, 2)

# prepare plot theme
plot_theme <- theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "gray60", size = 0.2),
    strip.text = element_text(size = 5),
    axis.text.x = element_text(size = 5), # , previous value: 6, angle = -45, hjust = -0.2),
    axis.text.y = element_text(size = 5), # previous value: 6,
    axis.line = element_line(colour = "black", size = 0.2),
    axis.ticks = element_line(size = 0.1), # previous value: 0.2
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 5), # previous value: 7
    legend.title = element_text(size = 6), # previous value: 8
    legend.position = "top",
    strip.background = element_rect(fill = "gray80", color = "black", size = 0.2),
    panel.background = element_rect(size = 0.2),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# create labs for plots
k.labs <- c("k = 5", "k = 10", "k = 30")
names(k.labs) <- c("5", "10", "30")

n.labs <- c("n = 50", "n = 100", "n = 200", "n = 500")
names(n.labs) <- c("50", "100", "200", "500")

or.labs <- c("logOR = -0.69", "logOR = 0", "logOR = 0.69")
names(or.labs) <- c("-0.69", "0", "0.69")

tau.labs <- c(paste(expression(tau), "= 0"), paste(expression(tau), "= 0.3"), 
              paste(expression(tau), "= 0.6"), paste(expression(tau), "= 1"))
names(tau.labs) <- c("0", "0.3", "0.6", "1")

p0.labs <- c("p0 = 0.01", "p0 = 0.05")
names(p0.labs) <- c("0.01", "0.05")

# restrict results_full to the models which are part of this manuscript:
results_full <- subset(results_full, model %in% c(
  "sj", "sj_knha", "mh", "binom_fi", "binom_fi2",
  "binom_ri", "binom_ri2", "betabinom", "nchg", "condbinom"
))

results_full$model <- factor(results_full$model,
  levels = c(
    "sj", "sj_knha", "mh", "binom_fi", "binom_fi2",
    "binom_ri", "binom_ri2", "betabinom", "nchg", "condbinom"
  ),
  labels = c(
    "IV_SJ", "IV_SJ_HKSJ", "MH", "BinFI", "BinFI2",
    "BinRI", "BinRI2", "Betabin", "Nchg", "Condbin"
  )
)

# DESCRIPTIVES ----------------------------------------------------------------------------

# Average percentage of single- and double-zero studies
# the following lines produce figures 1 - 3 of the online supplement:
# Figure 1 is named "avg00_k5.pdf"
# Figure 2 is named "avg00_k10.pdf"
# Figure 3 is named "avg00_k30.pdf"
for (k in c(5, 10, 30)) {
  plot_szdz <- results_full %>%
    filter(model == "IV_SJ", n_studies == k) %>%
    pivot_longer(
      cols = mean_n_dz:mean_n_sz_control,
      names_to = "group",
      values_to = "total",
      names_prefix = "mean_n_"
    ) %>%
    mutate(
      percentage = total / n_studies,
      group = factor(group,
        levels = c("dz", "sz_control", "sz_treat"),
        labels = c("double-zero", "single-zero (control)", "single-zero (treat)")
      )
    ) %>%
    ggplot(aes(x = factor(true_tau), y = percentage, group = group, 
               fill = group, colour = group)) +
    geom_bar(width = 0.3, position = position_dodge(width = 1), stat = "identity") +
    scale_y_continuous(position = "right") +
    facet_nested(true_logOR + dgm ~ p_control + n_sample,
      switch = "y",
      labeller = labeller(p_control = p0.labs, true_logOR = or.labs, 
                          n_sample = n.labs, n_studies = k.labs)
    ) +
    xlab("tau") +
    ylab("Average percentage of zero-studies") +
    guides(color = guide_legend("Double-zero/Single-zero: ")) +
    guides(fill = guide_legend("Double-zero/Single-zero: ")) +
    plot_theme

  plot_szdz

  plot_name <- paste0(plot_dir, "/avg00_k", k, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")
}


# Total number of simulation trials with zero-arm
# the following lines produce figures 4-6 of the online supplement:
# Figure 4 is named "zeroarm_k5.pdf"
# Figure 5 is named "zeroarm_k10.pdf"
# Figure 6 is named "zeroarm_k30.pdf"
for (k in c(5, 10, 30)) {
  plot_zeroarm <- results_full %>%
    filter(model == "IV_SJ", n_studies == k) %>%
    pivot_longer(
      cols = sum_all_dz:sum_all_sz_control,
      names_to = "group",
      values_to = "total",
      names_prefix = "sum_all_"
    ) %>%
    mutate(
      group = factor(group,
        levels = c("dz", "sz_control", "sz_treat"),
        labels = c("both", "control", "treat")
      )
    ) %>%
    ggplot(aes(x = factor(true_tau), y = total, group = group, 
               fill = group, colour = group)) +
    geom_bar(width = 0.3, position = position_dodge(width = 1), stat = "identity") +
    scale_y_continuous(position = "right", limits = c(0, 1000)) +
    facet_nested(true_logOR + dgm ~ p_control + n_sample,
      switch = "y",
      labeller = labeller(p_control = p0.labs, true_logOR = or.labs, 
                          n_sample = n.labs, n_studies = k.labs)
    ) +
    xlab("tau") +
    ylab("Number of simulation runs") +
    guides(color = guide_legend("Arm: ")) +
    guides(fill = guide_legend("Arm: ")) +
    plot_theme

  plot_zeroarm

  plot_name <- paste0(plot_dir, "/zeroarm_k", k, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")
}

# Average total number of events per arm
# the following lines produce figures 7-9 of the online supplement:
# Figure 7 is named "avgnevents_k5.pdf"
# Figure 8 is named "avgnevents_k10.pdf"
# Figure 9 is named "avgnevents_k30.pdf"
for (k in c(5, 10, 30)) {
  plot_tne <- results_full %>%
    filter(model == "IV_SJ", n_studies == k) %>%
    pivot_longer(
      cols = mean_total_treat:mean_total_control,
      names_to = "group",
      values_to = "mean_no_events",
      names_prefix = "mean_total_"
    ) %>%
    ggplot(aes(x = factor(true_tau), y = mean_no_events, group = group, 
               fill = group, colour = group)) +
    geom_bar(width = 0.6, position = position_dodge(width = 1), stat = "identity") +
    scale_y_continuous(position = "right") +
    facet_nested(true_logOR + dgm ~ p_control + n_sample,
      switch = "y",
      labeller = labeller(p_control = p0.labs, true_logOR = or.labs, 
                          n_sample = n.labs, n_studies = k.labs)
    ) +
    xlab("tau") +
    ylab("Mean total number of events") +
    guides(color = guide_legend("Group: ")) +
    guides(fill = guide_legend("Group: ")) +
    plot_theme

  plot_tne

  plot_name <- paste0(plot_dir, "/avgnevents_k", k, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")
}


# PERFORMANCE MEASURES ------------------------------------------------------------------

# create plot design data frame ---------------------------------------------------------

plot_design <- expand.grid(
  p_control = c(0.01, 0.05),
  group_ratio = c(1, 2),
  n_studies = c(5, 10, 30),
  n_sample = c("n_sample1", "n_sample2")
)


# the following lines produce figures 10 - 143 of the online supplement.
# the plots will be named according to the following naming scheme:
# "criterion_grA_pcB_kC_nD_nE.pdf"
# where "criterion" is replaced by the respective performance criterion
# (estimability, meanbias, medianbias, mse, coverage, ciwidth),
# A is replaced by the randomization ratio for which the results are depicted (1, 2),
# B is replaced by the event probability in the control group for which the results 
# are depicted (01 for p0 = 0.01, 05 for p0 = 0.05),
# C is replaced by the number of studies for which the results are depicted (5, 10, 30),
# D is replaced by the first sample size for which the results are depicted (50, 100),
# E is replaced by the second sample size for which the results are depicted (200, 500). 

for (i in 1:nrow(plot_design)) {
  gr <- plot_design$group_ratio[i]
  p0 <- plot_design$p_control[i]
  pc <- ifelse(p0 == 0.01, "01", "05")
  n1 <- ifelse(plot_design$n_sample[i] == "n_sample1", 100, 50)
  n2 <- ifelse(plot_design$n_sample[i] == "n_sample1", 500, 200)
  k <- plot_design$n_studies[i]

  # ESTIMABILITY -------------------------------------------------------------------------

  plot_est <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = estimability, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("Estimability") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_est


  plot_name <- paste0(plot_dir, "/estimability_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")


  # MEAN BIAS ---------------------------------------------------------------------------

  plot_meanbias <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = mean_bias_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    geom_hline(yintercept = 0, size = 0.15) +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("Mean Bias") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_meanbias


  plot_name <- paste0(plot_dir, "/meanbias_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")


  # MEDIAN BIAS -------------------------------------------------------------------------

  plot_medianbias <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = median_bias_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    geom_hline(yintercept = 0, size = 0.15) +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("Median Bias") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_medianbias

  plot_name <- paste0(plot_dir, "/medianbias_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")


  # MSE --------------------------------------------------------------------------------

  plot_mse <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = mse_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(0, 2), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("MSE") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_mse


  plot_name <- paste0(plot_dir, "/mse_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")


  # COVERAGE ---------------------------------------------------------------------------

  plot_coverage <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = coverage_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    geom_hline(yintercept = 0.95, size = 0.15, linetype = "dashed") +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(0.7, 1), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("Coverage") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_coverage


  plot_name <- paste0(plot_dir, "/coverage_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")


  # Median CI Width --------------------------------------------------------------------

  plot_ciwidth <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = median_ciwidth_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(0, 4), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
      switch = "y", axes = "all",
      labeller = labeller(
        n_sample = n.labs, p_control = p0.labs,
        true_logOR = or.labs, true_tau = tau.labs,
        n_studies = k.labs
      )
    ) +
    xlab("tau") +
    ylab("Median CI Width") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme

  plot_ciwidth


  plot_name <- paste0(plot_dir, "/ciwidth_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")
  
  # MCSE ----------------------------------------------------------------------------
  
  plot_mcse <- results_full %>%
    filter(
      group_ratio == gr, p_control == p0, n_studies == k, n_sample %in% c(n1, n2)
    ) %>%
    ggplot(aes(x = factor(true_tau), y = mcse_logOR, fill = model)) +
    geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
    scale_fill_manual(values = c(
      "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
      "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
    )) +
    scale_y_continuous(limits = c(0, 0.2), expand = c(0, 0), oob = rescale_none, 
                       position = "right") +
    scale_x_discrete(breaks = c("0", "0.3", "0.6", "1"), expand = c(0.15, 0.15)) +
    facet_nested(true_logOR + dgm ~ n_studies + n_sample,
                 switch = "y", axes = "all",
                 labeller = labeller(
                   n_sample = n.labs, p_control = p0.labs,
                   true_logOR = or.labs, true_tau = tau.labs,
                   n_studies = k.labs
                 )
    ) +
    xlab("tau") +
    ylab("MCSE") +
    guides(fill = guide_legend("Model: ")) +
    plot_theme
  
  plot_mcse
  
  
  plot_name <- paste0(plot_dir, "/mcse_gr", gr, "_pc", pc, "_k", k, "_n", 
                      n1, "_n", n2, ".pdf")
  ggsave(plot_name, device = "pdf", dpi = 300, height = 13, width = 18, units = "cm")
  
}
