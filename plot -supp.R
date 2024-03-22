##
## PRINT THE RESULTS
##

library(ggplot2)
library(gridExtra)
library(kableExtra)

svp <- F
svt <- F
## LOAD R DATA
load("Rdata/example-supp-bias.RData")

tab1_all[1,9] <- 1
## Pnmax = p
p1 <- ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(mapping=aes(ymin=BN.mu.lb,ymax=BN.mu.ub, colour = "BN-GLMM"), fill="red",  alpha=0.1, lty=2)+
  geom_ribbon(mapping=aes(ymin=NN.mu.lb,ymax=NN.mu.ub, colour = "NN model"),fill="grey50", alpha=0.1, lty=2)+
  geom_line(aes(y = BN.mu, colour="BN-GLMM"), lty=1, size=1) +
  geom_line(aes(y = NN.mu, colour="NN model"),lty=1, size=1) +
  geom_point(aes(y = BN.mu, colour="BN-GLMM"), size=2) +
  geom_point(aes(y = NN.mu, colour="NN model"), size=2) +
  scale_x_reverse(n.breaks = 10, name="Probability of publishing studies with smallest sample size / largest SE") +
  scale_y_continuous(limits = c(-7,0), name = "logit-transformed proportion") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.85), 
        legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Estimates with 95% CI",
                      breaks = c("BN-GLMM","NN model"),
                      values = c("red", "grey50"),
                      guide = guide_legend(
                        override.aes = 
                          list(lty = c(1,1), 
                               size = c(1,1),
                               fill = c("red", "grey50")
                               )))

if(svp) ggsave(filename = "plot-supp.eps", plot = p1, device = cairo_ps, width = 8, height = 8) 


## TABLE
load("Rdata/example-supp-bias.RData")

taba <- data.frame(
  p = sprintf("(%.1f, %.2f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  BN = sprintf("%.2f (%.2f, %.2f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.2f (%.2f)", c(lnOR_bn$tau[1], lnOR_copas_BNGLMM[5,]), c(lnOR_bn$tau[2], lnOR_copas_BNGLMM[6,])),
  BN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_BNGLMM[7,]), c(NA, lnOR_copas_BNGLMM[8,])),
  Mc = tab1_all$M.c,
  NN = sprintf("%.2f (%.2f, %.2f)", tab1_all$NN.mu, tab1_all$NN.mu.lb, tab1_all$NN.mu.ub),
  NN.tau = sprintf("%.2f (%.2f)", c(lnOR_nn$tau[1], lnOR_copas_NNLMM[5,]), c(lnOR_nn$tau[2], lnOR_copas_NNLMM[6,])),
  NN.rho = sprintf("%.2f (%.2f)", c(NA, lnOR_copas_NNLMM[7,]), c(NA, lnOR_copas_NNLMM[8,]))
)
taba[1,c(1,5,9)] <- ""

colnames(taba) <- c("$(\\Pmin, \\Pmax)$", "$M$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                    "$M$", "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")

if(svt) sink("tab-supp.tex")
kbl(taba, 
    format = ifelse(svt, "latex", "html"),
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab")%>% 
  add_header_above(c("", "Proposed method (BN-GLMM)" = 4, "Copas-Heckman-type selection model" = 4))%>%
  footnote(general = "$M$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error;
           NaN indicates ``not a value'', an undefined value.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
if(svt) sink()