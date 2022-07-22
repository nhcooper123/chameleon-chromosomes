# Make phylogeny figure(s)

# Load libraries
library(tidyverse)
library(here)
library(phytools)
library(treeplyr)
source(here("chromEvol/ChromEvol.R")) #install this from source
library(ggtree)

# Read in the results for the best fitting models
results <- read.ce(here("chromEvol/results/02D_18-root-exclude_out/CONST_RATE_NO_DUPL"))
results_all <- read.ce(here("chromEvol/results/03D_estimate-exclude_out/CONST_RATE_DEMI_EST"))
results_est <- read.ce(here("chromEvol/results/03D_estimate-exclude_out/CONST_RATE_NO_DUPL"))

# Read in haploid number data
taxa <- read_csv(here("data/chromosome-data-species.csv"))
taxa <- as.data.frame(taxa)
taxa <- filter(taxa, tips != "Riebrev" & tips != "Leiolep")

# Read in the tree
mytree <- read.tree(here("chromEvol/tree_for_chromevol_exclude"))

#---------------------------------------------------------
# Prep for plotting
#-----------------------------------------------------------
# Extract required data for plotting from ChromEvol object
# pp gives the number with the highest posterior probability
taxa_anc <- results[[1]][[1]]
taxa_tip <- results[[1]][[2]]
taxa_pp <- results[[1]][[3]]

# Create list of colours
mycolours <- c("#A4C61A","#208EA3", "#FD817D", "#4178BC", "#7A71F6",
                "#EECC16", "#E8384F", "#E37CFF", "#37A862", "#EA4E9D")


#"#E8384F" = RED
#"#208EA3" = AQUA
#"#FD817D" = ORANGE
#"#4178BC" = BLUE
#"#FDAE33" = YELLOW ORANGE
#"#7A71F6" = INDIGO
#"#EECC16" = YELLOW
#"#AA71FF" = PURPLE
#"#A4C61A" = YELLOW GREEN
#"#E37CFF" = MAGENTA
#"#62BB35" = GREEN
#"#EA4E9D" = HOT PINK
#"#37A862" = BLUE GREEN

plot(mytree, type = "fan",
     show.tip.label = FALSE, no.margin = TRUE)
#tiplabels(pch = 16, col = tip_col, offset = 0.5)
tiplabels(pie = taxa_tip, piecol = mycolours, offset = 5, cex = 0.6)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = taxa_anc, piecol = mycolours, cex = 0.8)
nodelabels(text = taxa_pp, frame = "circle", bg = "white", cex = 0.6)

# Save the plot using the export function once it is the correct size.

#--------------------------------------------------
# ALL MODELS
#--------------------------------------------------

all_anc <- results_all[[1]][[1]]
all_tip <- results_all[[1]][[2]]
all_pp <- results_all[[1]][[3]]
all_phy <- results_all[[3]]

plot(mytree, type = "fan",
     show.tip.label = FALSE, no.margin = TRUE)
#tiplabels(pch = 16, col = tip_col, offset = 0.5)
tiplabels(pie = all_tip, piecol = mycolours, offset = 5, cex = 0.6)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = all_anc, piecol = mycolours, cex = 0.8)
nodelabels(text = all_pp, frame = "circle", bg = "white", cex = 0.6)

# Save the plot using the export function once it is the correct size.

#--------------------------------------------------
# ESTIMATED
#--------------------------------------------------

est_anc <- results_est[[1]][[1]]
est_tip <- results_est[[1]][[2]]
est_pp <- results_est[[1]][[3]]
est_phy <- results_est[[3]]

plot(mytree, type = "fan",
     show.tip.label = FALSE, no.margin = TRUE)
#tiplabels(pch = 16, col = tip_col, offset = 0.5)
tiplabels(pie = est_tip, piecol = mycolours, offset = 5, cex = 0.6)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = est_anc, piecol = mycolours, cex = 0.8)
nodelabels(text = est_pp, frame = "circle", bg = "white", cex = 0.6)

# Save the plot using the export function once it is the correct size.

