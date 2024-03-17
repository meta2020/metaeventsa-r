##
## PRINT THE RESULTS
##

library(ggplot2)
library(gridExtra)
library(kableExtra)

svp <- T
svt <- T
## LOAD R DATA
load("Rdata/example-bias1.RData")

tab1_all[1,12] <- 1
## Pnmax = p
p1 <- ggplot(tab1_all, aes(x = pnmax)) +
  geom_ribbon(mapping=aes(ymin=HN.mu.lb,ymax=HN.mu.ub, colour = "HN-GLMM"), fill="red", alpha=0.1, lty=2)+
  geom_ribbon(mapping=aes(ymin=BN.mu.lb,ymax=BN.mu.ub, colour = "BN-GLMM"), fill="steelblue",  alpha=0.1, lty=2)+
  geom_ribbon(mapping=aes(ymin=NN.mu.lb,ymax=NN.mu.ub, colour = "NN model"),fill="grey50", alpha=0.1, lty=2)+
  geom_line(aes(y = HN.mu, colour="HN-GLMM"), lty=1, size=1) +
  geom_line(aes(y = BN.mu, colour="BN-GLMM"), lty=1, size=1) +
  geom_line(aes(y = NN.mu, colour="NN model"),lty=1, size=1) +
  geom_point(aes(y = HN.mu, colour="HN-GLMM"), size=2) +
  geom_point(aes(y = BN.mu, colour="BN-GLMM"), size=2) +
  geom_point(aes(y = NN.mu, colour="NN model"), size=2) +
  scale_x_reverse(n.breaks = 10, name="Probability of publishing studies with largerst sample size / smallest SE") +
  scale_y_continuous(limits = c(-3,0), name = "lnOR") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.15), 
        legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Estimates with 95% CI",
                      breaks = c("HN-GLMM", "BN-GLMM","NN model"),
                      values = c("red", "steelblue", "grey50"),
                      guide = guide_legend(
                        override.aes = 
                          list(lty = c(1,1,1), 
                               size = c(1,1,1),
                               fill = c("red", "steelblue", "grey50")
                               )))+
  ggtitle("B")


## LOAD R DATA
load("Rdata/example-bias2.RData")

tab2_all[1,12] <- 1
## Pnmax = p
p2 <- ggplot(tab2_all, aes(x = pnmin)) +
  geom_ribbon(mapping=aes(ymin=HN.mu.lb,ymax=HN.mu.ub, colour = "HN-GLMM"), fill="red", alpha=0.1, lty=2)+
  geom_ribbon(mapping=aes(ymin=BN.mu.lb,ymax=BN.mu.ub, colour = "BN-GLMM"), fill="steelblue",  alpha=0.1, lty=2)+
  geom_ribbon(mapping=aes(ymin=NN.mu.lb,ymax=NN.mu.ub, colour = "NN model"),fill="grey50", alpha=0.1, lty=2)+
  geom_line(aes(y = HN.mu, colour="HN-GLMM"), lty=1, size=1) +
  geom_line(aes(y = BN.mu, colour="BN-GLMM"), lty=1, size=1) +
  geom_line(aes(y = NN.mu, colour="NN model"),lty=1, size=1) +
  geom_point(aes(y = HN.mu, colour="HN-GLMM"), size=2) +
  geom_point(aes(y = BN.mu, colour="BN-GLMM"), size=2) +
  geom_point(aes(y = NN.mu, colour="NN model"), size=2) +
  scale_x_reverse(n.breaks = 10, name="Probability of publishing studies with smallest sample size / largest SE") +
  scale_y_continuous(limits = c(-3,0), name = "lnOR") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.15), 
        legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Estimates with 95% CI",
                      breaks = c("HN-GLMM", "BN-GLMM","NN model"),
                      values = c("red", "steelblue", "grey50"),
                      guide = guide_legend(
                        override.aes = 
                          list(lty = c(1,1,1), 
                               size = c(1,1,1),
                               fill = c("red", "steelblue", "grey50")
                          )))+
  ggtitle("A")

plot <- grid.arrange(p2, p1, ncol=2)
if(svp) ggsave(filename = "plot.eps", plot = plot, device = cairo_ps, width = 12, height = 6) 


## TABLE
load("Rdata/example-bias1.RData")

taba <- data.frame(
  p = sprintf("(%.1f, %.2f)", tab1_all$pnmax, tab1_all$pnmin),
  Mp = tab1_all$M.p,
  HN = sprintf("%.2f (%.2f, %.2f)", tab1_all$HN.mu, tab1_all$HN.mu.lb, tab1_all$HN.mu.ub),
  HN.tau = sprintf("%.2f (%.2f)", c(lnOR_hn$tau[1], lnOR_copas_HNGLMM[5,]), c(lnOR_hn$tau[2], lnOR_copas_HNGLMM[6,])),
  HN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_HNGLMM[7,]), c(NA, lnOR_copas_HNGLMM[8,])),
  BN = sprintf("%.2f (%.2f, %.2f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.2f (%.2f)", c(lnOR_bn$tau[1], lnOR_copas_BNGLMM[5,]), c(lnOR_bn$tau[2], lnOR_copas_BNGLMM[6,])),
  BN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_BNGLMM[7,]), c(NA, lnOR_copas_BNGLMM[8,])),
  Mc = tab1_all$M.c,
  NN = sprintf("%.2f (%.2f, %.2f)", tab1_all$NN.mu, tab1_all$NN.mu.lb, tab1_all$NN.mu.ub),
  NN.tau = sprintf("%.2f (%.2f)", c(lnOR_nn$tau[1], lnOR_copas_NNLMM[5,]), c(lnOR_nn$tau[2], lnOR_copas_NNLMM[6,])),
  NN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_NNLMM[7,]), c(NA, lnOR_copas_NNLMM[8,]))
)
taba[1,c(1,5,8,12)] <- ""

load("Rdata/example-bias2.RData")

tabb <- data.frame(
  p = sprintf("(%.2f, %.1f)", tab2_all$pnmax, tab2_all$pnmin),
  Mp = tab2_all$M.p,
  HN = sprintf("%.2f (%.2f, %.2f)", tab2_all$HN.mu, tab2_all$HN.mu.lb, tab2_all$HN.mu.ub),
  HN.tau = sprintf("%.2f (%.2f)", c(lnOR_hn$tau[1], lnOR_copas_HNGLMM[5,]), c(lnOR_hn$tau[2], lnOR_copas_HNGLMM[6,])),
  HN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_HNGLMM[7,]), c(NA, lnOR_copas_HNGLMM[8,])),
  BN = sprintf("%.2f (%.2f, %.2f)", tab2_all$BN.mu, tab2_all$BN.mu.lb, tab2_all$BN.mu.ub),
  BN.tau = sprintf("%.2f (%.2f)", c(lnOR_bn$tau[1], lnOR_copas_BNGLMM[5,]), c(lnOR_bn$tau[2], lnOR_copas_BNGLMM[6,])),
  BN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_BNGLMM[7,]), c(NA, lnOR_copas_BNGLMM[8,])),
  Mc = tab2_all$M.c,
  NN = sprintf("%.2f (%.2f, %.2f)", tab2_all$NN.mu, tab2_all$NN.mu.lb, tab2_all$NN.mu.ub),
  NN.tau = sprintf("%.2f (%.2f)", c(lnOR_nn$tau[1], lnOR_copas_NNLMM[5,]), c(lnOR_nn$tau[2], lnOR_copas_NNLMM[6,])),
  NN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_NNLMM[7,]), c(NA, lnOR_copas_NNLMM[8,]))
)
tabb[1,c(1,5,8,12)] <- ""

tab_all <- rbind(tabb, taba)
colnames(tab_all) <- c("$(\\Pmin, \\Pmax)$", "$M$", 
                       "$\\theta$ (95% CI)", "$\\tau (SE)$", "$\\rho (SE)$", 
                       "$\\theta$ (95% CI)", "$\\tau (SE)$", "$\\rho (SE)$", 
                       "$M$", "$\\theta$ (95% CI)", "$\\tau (SE)$", "$\\rho (SE)$")

if(svt) sink("tab.tex")
kbl(tab_all, 
    format = ifelse(svt, "latex", "html"),
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab")%>% 
  add_header_above(c(" ", "Proposed method (HN-GLMM)" = 4, "Proposed method (BN-GLMM)" = 3, "Copas-Heckman-type selection model" = 4))%>%
  footnote(general = "$M$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error;
           NaN indicates ``not a value'', an undefined value.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
if(svt) sink()