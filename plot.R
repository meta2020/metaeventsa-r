##
## PRINT THE RESULTS
##
rm(list=ls())
library(ggplot2)
library(gridExtra)
library(kableExtra)


## LOAD R DATA
load("Rdata/example-bias.RData")

## Pnmax = p
p1 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = HN.mu.lb, ymax = HN.mu.ub), alpha = 0.1, fill = "#d7191c", na.rm = TRUE) + 
  geom_ribbon(aes(ymin = BN.mu.lb, ymax = BN.mu.ub), alpha = 0.1, fill = "#fdae61", na.rm = TRUE) + 
  geom_ribbon(aes(ymin = CH.mu.lb, ymax = CH.mu.ub), alpha = 0.1, fill = "#2c7bb6", na.rm = TRUE) + 
  geom_ribbon(aes(ymin = CN.mu.lb, ymax = CN.mu.ub), alpha = 0.1, fill = "#abd9e9", na.rm = TRUE) + 
  geom_line(aes(y = HN.mu.lb, colour="HN-GLMM"), lty=2, size=1) +
  geom_line(aes(y = HN.mu.ub, colour="HN-GLMM"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.lb, colour="BN-GLMM"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.ub, colour="BN-GLMM"), lty=2, size=1) +
  geom_line(aes(y = CH.mu.lb, colour="CH model"),lty=2, size=1) +
  geom_line(aes(y = CH.mu.ub, colour="CH model"),lty=2, size=1) +
  geom_line(aes(y = CN.mu.lb, colour="CN model"),lty=2, size=1) +
  geom_line(aes(y = CN.mu.ub, colour="CN model"),lty=2, size=1) +
  geom_point(aes(y = HN.mu, colour="HN-GLMM"), size=3) +
  geom_point(aes(y = BN.mu, colour="BN-GLMM"), size=3) +
  geom_point(aes(y = CH.mu, colour="CH model"), size=3) +
  geom_point(aes(y = CN.mu, colour="CN model"), size=3) +
  geom_line(aes(y = HN.mu, colour="HN-GLMM"), lty=1, size=1.5) +
  geom_line(aes(y = BN.mu, colour="BN-GLMM"), lty=1, size=1.5) +
  geom_line(aes(y = CH.mu, colour="CH model"),lty=1, size=1.5) +
  geom_line(aes(y = CN.mu, colour="CN model"),lty=1, size=1.5) +
  scale_x_reverse(n.breaks = 10, name="Probability of publishing studies with smallest sample size / largest SE") +
  scale_y_continuous(limits = c(-3,0), name = "lnOR", n.breaks = 10) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        # legend.position = "inside",
        legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Estimates with 95% CI",
                      breaks = c("HN-GLMM", "BN-GLMM","CH model","CN model"),
                      values = c("#d7191c", "#fdae61","#2c7bb6","#abd9e9"),
                      guide = guide_legend(
                        override.aes = 
                          list(lty = c(1,1,1,1), 
                               size = c(1,1,1,1)
                               )
                        ))

ggsave(filename = "plot.eps", plot = p1, device = cairo_ps, width = 8, height = 5) 


## TABLE
load("Rdata/example-bias.RData")

taba = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  HN = sprintf("%.3f (%.3f, %.3f)", tab1_all$HN.mu, tab1_all$HN.mu.lb, tab1_all$HN.mu.ub),
  HN.tau = sprintf("%.3f (%.3f)", lnOR_COPAS_HNGLMM[5,], lnOR_COPAS_HNGLMM[6,]),
  HN.rho = sprintf("%.3f (%.3f)", lnOR_COPAS_HNGLMM[7,], lnOR_COPAS_HNGLMM[8,]),
  BN = sprintf("%.3f (%.3f, %.3f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.3f (%.3f)", lnOR_COPAS_BNGLMM[5,], lnOR_COPAS_BNGLMM[6,]),
  BN.rho = sprintf("%.3f (%.3f)", lnOR_COPAS_BNGLMM[7,], lnOR_COPAS_BNGLMM[8,]),
  CN = sprintf("%.3f (%.3f, %.3f)", tab1_all$CN.mu, tab1_all$CN.mu.lb, tab1_all$CN.mu.ub),
  CN.tau = sprintf("%.3f (%.3f)", lnOR_COPAS1999[5,], lnOR_COPAS1999[6,]),
  CN.rho = sprintf("%.3f (%.3f)", lnOR_COPAS1999[7,], lnOR_COPAS1999[8,]),
  Mc = tab1_all$M.c,
  CH = sprintf("%.3f (%.3f, %.3f)", tab1_all$CH.mu, tab1_all$CH.mu.lb, tab1_all$CH.mu.ub),
  CH.tau = sprintf("%.3f (%.3f)", lnOR_COPAS2000[5,], lnOR_COPAS2000[6,]),
  CH.rho = sprintf("%.3f (%.3f)", lnOR_COPAS2000[7,], lnOR_COPAS2000[8,])
)

colnames(taba) = c("$(\\Pmin, \\Pmax)$", "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                   "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")
taba1 = taba[,1:8]
colnames(taba1) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                   "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")
taba2 = taba[,c(1,2,9:15)]
colnames(taba2) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                    "$M2$",
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")

kbl(taba1, 
    format = "html",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab")%>% 
  add_header_above(c(" ","", 
                     "Proposed method (HN-GLMM)" = 3, 
                     "Proposed method (BN-GLMM)" = 3
                     ))%>%
  footnote(general = "$M1$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

kbl(taba2, 
    format = "html",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab")%>% 
  add_header_above(c(" ","", 
                     "Copas-N selection model" = 3,"",
                     "Copas-Heckman selection model" = 3))%>%
  footnote(general = "$M1$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
