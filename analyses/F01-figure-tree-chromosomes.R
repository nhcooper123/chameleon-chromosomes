# Make phylogeny figure(s)

# Load libraries
library(tidyverse)
library(here)
library(phytools)
library(treeplyr)
library(ChromEvol) #install this from source
library(here)
library(ggtree)

# Read in the results
results_taxa <- read.ce(here("chromEvol/05_chromEvol_out_18-root_exclude/CONST_RATE_NO_DUPL"))
results_species <- read.ce(here("chromEvol/11_chromEvol_out_18-root_exclude_species/CONST_RATE_NO_DUPL"))

# Read in haploid number data
taxa <- read_csv(here("data/chromosome-data.csv"))
taxa <- as.data.frame(taxa)
taxa <- filter(taxa, tips != "Riebrev" & tips != "Leiolep")

# Read in the tree
mytree <- read.tree(here("chromEvol/tree_for_chromevol_exclude"))

#---------------------------------------------------------
# Prep for plotting
#-----------------------------------------------------------
# Extract required data for plotting from ChromEvol object
# pp gives the number with the highest posterior probability
taxa_anc <- results_taxa[[1]][[1]]
taxa_tip <- results_taxa[[1]][[2]]
taxa_pp <- results_taxa[[1]][[3]]

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
# SPECIES LEVEL
#--------------------------------------------------

# Read in the tree
mytree <- read.tree(here("chromEvol/tree_for_chromevol_exclude_species"))

# Remove excess for species data
species <- 
  taxa %>% 
  mutate(dup = duplicated(taxa$Binomial_Marcello)) %>%
  filter(dup != TRUE)

species_anc <- results_species[[1]][[1]]
species_tip <- results_species[[1]][[2]]
species_pp <- results_species[[1]][[3]]
species_phy <- results_species[[3]]

plot(mytree, type = "fan",
     show.tip.label = FALSE, no.margin = TRUE)
#tiplabels(pch = 16, col = tip_col, offset = 0.5)
tiplabels(pie = species_tip, piecol = mycolours, offset = 5, cex = 0.6)
tiplabels(text = species$haploidn, frame = "none", cex = 0.6, 
          offset = 5)
nodelabels(pie = species_anc, piecol = mycolours, cex = 0.8)
nodelabels(text = species_pp, frame = "circle", bg = "white", cex = 0.6)

# Save the plot using the export function once it is the correct size.
