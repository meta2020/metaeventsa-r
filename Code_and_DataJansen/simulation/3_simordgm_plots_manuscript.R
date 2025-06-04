#---------------------------------------------------------------------------------------------
# PLOTS (MANUSCRIPT)
#---------------------------------------------------------------------------------------------

# this script produces the plots for the manuscript.
# it requires that the file "simordgm_summary.RData", which contains the summarized
# simulation results, exists in the results folder in the working directory.
# this script creates a directory called plots_manuscript in which the
# figures are saved as eps files.

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
dir.create("plots_manuscript")
plot_dir <- "./plots_manuscript"

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

tau.labs <- c(paste(expression(tau), "= 0"), paste(expression(tau), "= 0.3"), paste(expression(tau), "= 0.6"), paste(expression(tau), "= 1"))
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

# ESTIMABILITY -----------------------------------------------------------------------------------

# Figure 1

figure1 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 5, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = estimability, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure1

plot_name <- paste0(plot_dir, "/figure1.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")

# Figure 2

figure2 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 30, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = estimability, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure2

plot_name <- paste0(plot_dir, "/figure2.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")




# MEDIAN BIAS ------------------------------------------------------------------------------------

# Figure 3

figure3 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 5, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = median_bias_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  geom_hline(yintercept = 0, size = 0.15) +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure3

plot_name <- paste0(plot_dir, "/figure3.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")


# Figure 4

figure4 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 30, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = median_bias_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  geom_hline(yintercept = 0, size = 0.15) +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure4

plot_name <- paste0(plot_dir, "/figure4.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")

# MSE --------------------------------------------------------------------------------------------

# Figure 5

figure5 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 5, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = mse_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure5


plot_name <- paste0(plot_dir, "/figure5.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")

# Figure 6

figure6 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 30, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = mse_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure6


plot_name <- paste0(plot_dir, "/figure6.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")


# COVERAGE --------------------------------------------------------------------------------------

# Figure 7

figure7 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 5, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = coverage_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  geom_hline(yintercept = 0.95, size = 0.15, linetype = "dashed") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0.7, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure7

plot_name <- paste0(plot_dir, "/figure7.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")


# Figure 8

figure8 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 30, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = coverage_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  geom_hline(yintercept = 0.95, size = 0.15, linetype = "dashed") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0.7, 1), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure8

plot_name <- paste0(plot_dir, "/figure8.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")

# Median CI Width --------------------------------------------------------------------------------------------

# Figure 9

figure9 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 5, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = median_ciwidth_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure9


plot_name <- paste0(plot_dir, "/figure9.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")


# Figure 10

figure10 <- results_full %>%
  filter(
    group_ratio == 1, p_control == 0.01, n_studies == 30, n_sample %in% c(100, 500)
  ) %>%
  ggplot(aes(x = factor(true_tau), y = median_ciwidth_logOR, fill = model)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "identity") +
  scale_fill_manual(values = c(
    "darkblue", "cadetblue2", "grey", "darkgreen", "darkolivegreen2",
    "darkred", "coral2", "darkorange", "purple4", "mediumpurple2"
  )) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0), oob = rescale_none, position = "right") +
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

figure10


plot_name <- paste0(plot_dir, "/figure10.eps")
ggsave(plot_name, device = "eps", dpi = 300, height = 13, width = 18, units = "cm")

