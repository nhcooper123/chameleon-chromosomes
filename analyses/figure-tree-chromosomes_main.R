# Make phylogeny figure(s)

# Load libraries
library(tidyverse)
library(here)
library(phytools)
source(here("chromEvol/ChromEvol.R")) #install this from source
library(ggtree)

# Read in the results for the best fitting models
results <- read.ce(here("chromEvol/results/02D_18-root-exclude_out/CONST_RATE_NO_DUPL"))

# Read in haploid number data
taxa <- read_csv(here("data/chromosome-data-species.csv"))
taxa <- filter(taxa, tips != "Riebrev" & tips != "Leiolep")

# Read in the tree
mytree <- read.tree(here("chromEvol/tree_for_chromevol_exclude"))

# Reorder data so it is the same order as the tree tip labels
taxa <- taxa[match(mytree$tip.label, taxa$tips), ]
# Replace tips with Binomials
mytree$tip.label <- taxa$Binomial
#---------------------------------------------------------
# Prep for plotting
#-----------------------------------------------------------
# Extract required data for plotting from ChromEvol object
# pp gives the number with the highest posterior probability
taxa_anc <- results[[1]][[1]]
taxa_tip <- results[[1]][[2]]
taxa_pp <- results[[1]][[3]]

# Create list of colours (colourblind friendly from: https://personal.sron.nl/~pault/#sec:qualitative)
mycolours <- c("#332288", "#33BBEE", "#44AA99", "#117733",
               "#62BB35", "#FDAE33", "#CC6677", "#EE3377", 
               "#AA4499", "#AA71FF")

taxa <- 
  taxa %>%
  mutate(colours = case_when(haploidn == 10 ~ "#332288",
                             haploidn == 11 ~ "#33BBEE",
                             haploidn == 12 ~ "#44AA99",
                             haploidn == 13 ~ "#117733",
                             haploidn == 14 ~ "#62BB35",
                             haploidn == 15 ~ "#FDAE33",
                             haploidn == 16 ~ "#CC6677",
                             haploidn == 17 ~ "#EE3377",
                             haploidn == 18 ~ "#AA4499",
                             haploidn == 19 ~ "#AA71FF"))

taxa_pp <- 
  data.frame(taxa_pp) %>%
  mutate(colours = case_when(taxa_pp == 10 ~ "#332288",
                             taxa_pp == 11 ~ "#33BBEE",
                             taxa_pp == 12 ~ "#44AA99",
                             taxa_pp == 13 ~ "#117733",
                             taxa_pp == 14 ~ "#62BB35",
                             taxa_pp == 15 ~ "#FDAE33",
                             taxa_pp == 16 ~ "#CC6677",
                             taxa_pp == 17 ~ "#EE3377",
                             taxa_pp == 18 ~ "#AA4499",
                             taxa_pp == 19 ~ "#AA71FF"))
#-----------------------------------------------------------
# plot (minus pies)
plot(mytree, 
     show.tip.label = TRUE, no.margin = TRUE, 
     tip.color = taxa$colours, cex = 0.5, font = 3, 
     label.offset = 5)
tiplabels(pch = 16, col = taxa$colours, offset = 0.5)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.5, 
         offset = 3)
#nodelabels(pie = taxa_anc, piecol = mycolours, cex = 0.6)
nodelabels(pch = 16, col = taxa_pp$colours, cex = 3)
nodelabels(text = taxa_pp$taxa_pp, frame = "circle", bg = "white", cex = 0.5)

# Save the plot using the export function once it is the correct size.
# width 600/700
#-----------------------------------------------------------
# No colours for tips
plot(mytree, 
     show.tip.label = TRUE, no.margin = TRUE, 
     tip.color = "black", cex = 0.5, font = 3, 
     label.offset = 5)
tiplabels(pch = 16, col = taxa$colours, offset = 0.5)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.5, 
          offset = 3)
#nodelabels(pie = taxa_anc, piecol = mycolours, cex = 0.6)
nodelabels(pch = 16, col = taxa_pp$colours, cex = 3)
nodelabels(text = taxa_pp$taxa_pp, frame = "circle", bg = "white", cex = 0.5)
#-----------------------------------------------------------
# FAN
par(mar = c(0.01, .01, .01, .01))
plot(mytree, 
     show.tip.label = TRUE, no.margin = FALSE, 
     tip.color = taxa$colours, cex = 0.7, font = 3, 
     label.offset = 5, type = "fan")
tiplabels(pch = 16, col = taxa$colours, offset = 0.5)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.5, 
          offset = 3)
#nodelabels(pie = taxa_anc, piecol = mycolours, cex = 0.6)
nodelabels(pch = 16, col = taxa_pp$colours, cex = 4)
nodelabels(text = taxa_pp$taxa_pp, frame = "circle", bg = "white", cex = 0.7)
# Width = 950
#-----------------------------------------------------------
# With pies

# Numbers 1-9 grey, 10-19 colours, 20-29 - grey.
mycolours_pies <- c(rep("#DDDDDD",9),
                "#332288", "#33BBEE", "#44AA99", "#117733",
               "#62BB35", "#FDAE33", "#CC6677", "#EE3377", 
               "#AA4499", "#AA71FF", rep("#DDDDDD", 10))

# Fan
plot(mytree, type = "fan",
     show.tip.label = TRUE, no.margin = FALSE, cex = 0.7, font = 3, 
     label.offset = 10, tip.color = taxa$colours)
tiplabels(pie = taxa_tip, piecol = mycolours_pies, offset = 5, cex = 0.6)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = taxa_anc, piecol = mycolours_pies, cex = 0.8)
nodelabels(text = taxa_pp$taxa_pp, frame = "circle", bg = "white", cex = 0.6)

# Save the plot using the export function once it is the correct size.
# THIS IS MY FAVOURITE PLOT
plot(mytree, type = "fan",
     show.tip.label = TRUE, no.margin = FALSE, cex = 0.7, font = 3, 
     label.offset = 10)
tiplabels(pie = taxa_tip, piecol = mycolours_pies, offset = 5, cex = 0.6)
tiplabels(text = taxa$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = taxa_anc, piecol = mycolours_pies, cex = 0.8)
nodelabels(text = taxa_pp$taxa_pp, frame = "circle", bg = "white", cex = 0.6)