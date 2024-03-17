## [Reproducible R Codes] 

# Copas-Heckman-type sensitivity analysis for publication bias in rare-event meta-analysis under the framework of generalized linear mixed model


This folder contains reproducible R codes of simulation studies and re-analysis of the example data.

The following packages are used

- `cubature`, `MCMCpack`, `numDeriv` in the likelihoods

- `ggplot2`, `kableExtra`, `gridExtra` are used for plots and tables

If they are not installed, please install from R CRAN `install.packages("package_name")`.
 

## R codes for the estimation of the examples

- [Example in the main text](example.R)

- [Example in the supplementary file](example-supp.R)

## R codes for generating plots and tables

- [Example in the main text](plot.R)

- [Example in the supplementary file](plot-supp.R)

R codes of the models with/without publication bias are in folder "Rfunc"; 
the estimated results are saved in R data in folder "Rdata".