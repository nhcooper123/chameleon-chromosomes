# Get chromEvol outputs
# Load libraries
library(tidyverse)
library(here)
library(phytools)
source(here("chromEvol/ChromEvol.R")) #install this from source
library(ggtree)

#-----------------------------------------------
# Read in the results for the estimating models
#-----------------------------------------------
results1 <- read.ce(here("chromEvol/results/03A_estimate_noO_out/CONST_RATE_NO_DUPL"))
results2 <- read.ce(here("chromEvol/results/03B_estimate-exclude_out/CONST_RATE_NO_DUPL"))

results3 <- read.ce(here("chromEvol/OTUs/results/03A_estimate_noO_out/CONST_RATE_NO_DUPL"))
results4 <- read.ce(here("chromEvol/OTUs/results/03B_estimate-exclude_out/CONST_RATE_NO_DUPL"))

#-----------------------------------------------
# What is best root estimate in terms of pp?
#-----------------------------------------------
# pp gives the number with the highest posterior probability
results1[[1]][[3]][1]
results2[[1]][[3]][1]
results3[[1]][[3]][1]
results4[[1]][[3]][1]
