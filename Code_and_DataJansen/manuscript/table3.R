# Table 3, manuscript:
# Random-effects meta-analysis models for the odds ratio in the case of rare events
# under different data generating models: A simulation study

# this script reproduces Table 3 from the abovementioned manuscript.

# if it is desired to produce the latex table, the xtable package can be used:
library(xtable)

# specify conditions for which intervals are calculated:
conditions <- expand.grid(
  tau = c(0, 0.3, 0.6, 1),
  p0 = c(0.01, 0.05)
)

# calculate logit(p0)
conditions$alpha <- log(conditions$p0 / (1 - conditions$p0))

# calculate lower and upper bounds of the 95 % interval for the logits:
conditions$alpha.lb <- conditions$alpha - qnorm(0.975) * conditions$tau / sqrt(2)
conditions$alpha.ub <- conditions$alpha + qnorm(0.975) * conditions$tau / sqrt(2)

# backtransform lower and upper bound of the 95 % interval to probability scale:
conditions$p0.lb <- 1 / (1 + exp(-conditions$alpha.lb))
conditions$p0.ub <- 1 / (1 + exp(-conditions$alpha.ub))

# create character vector which contains both interval bounds:
conditions$interval <- paste0(round(conditions$p0.lb, 3), " - ", round(conditions$p0.ub, 3))

# create Table 3:
table3 <- rbind(
  conditions[conditions$p0 == 0.01 & conditions$tau > 0, ]$interval,
  conditions[conditions$p0 == 0.05 & conditions$tau > 0, ]$interval
)
rownames(table3) <- c("p_0 = 0.01", "p_0 = 0.05")
colnames(table3) <- c("\\tau = 0.3", "\\tau = 0.6", "\\tau = 1")

# show Table 3:
table3

# produce and print latex table:
xtable3 <- xtable(table3)
print(xtable3, sanitize.text.function = function(x) {
  x
})
