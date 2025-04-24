# Table B1 and Table B2,
# manuscript: Random-effects meta-analysis models for the odds ratio 
# in the case of rare events under different data generating models: A simulation study

# this script reproduces Table B1 and Table B2 from the abovementioned manuscript.

library(xtable)

setwd("~/example")

# Prepare data -------------------------------------------------------------------------------
load("data_thomas.RData")

# extract data for the outcomes suicidal ideation and aggression:
tableB1 <- subset(data_thomas, outcome == "Suicidal Ideation")
tableB2 <- subset(data_thomas, outcome == "Aggression")

# add information on study names to the tables:
studiesB1 <- c("Jorenby 2006", "Hughes 2011", "Bolliger 2011",
               "Tashkin 2011", "Williams 2012", "Rennard 2012",
               "Mitchell 2012", "Stein 2013", "Cinciripini 2013",
               "Meszaros 2013", "Anthenelli 2013", "Evins 2014",
               "Chengappa 2014", "Zhao 2011", "Brandon 2011",
               "Ebbert 2011", "Steinberg 2011", "Wong 2012", 
               "McClure 2013","Gonzales 2014")

studiesB2 <- c("Tonstad 2006", "Jorenby 2006", "Gonzales 2006",
               "Williams 2007", "Tsai 2007", "Faessel 2009",
               "Rigotti 2010", "Hughes 2011", "Garza 2011",
               "Bolliger 2011", "Tashkin 2011", "Williams 2012",
               "Rennard 2012", "Mitchell 2012", "Stein 2013",
               "Litten 2013", "Meszaros 2013", "Anthenelli 2013",
               "Chengappa 2014", "Gonzales 2014", "Brandon 2011",
               "Steinberg 2011"
               )

tableB1 <- cbind(studiesB1, tableB1[,c("event_treat", "n_treat", "event_control", "n_control")])
tableB2 <- cbind(studiesB2, tableB2[,c("event_treat", "n_treat", "event_control", "n_control")])

# make latex tables:
names(tableB1) <- names(tableB2) <- c("Study", "y_{i1}", "n_{i1}", "y_{i0}", "n_{i0}")

xtableB1 <- xtable(tableB1, digits = 0)
xtableB2 <- xtable(tableB2, digits = 0)

print(xtableB1, sanitize.text.function = function(x){x}, include.rownames = FALSE)
print(xtableB2, sanitize.text.function = function(x){x}, include.rownames = FALSE)
