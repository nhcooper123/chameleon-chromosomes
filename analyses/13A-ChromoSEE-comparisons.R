# Compare ancestral states on tree and dates

#----------------------------------
# Load libraries
#----------------------------------
library(RevGadgets)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ape)

#----------------------------------
# Get correct tip labels
#----------------------------------
# Read in haploid number data
mydata <- read_csv("data/chromosome-data-species1.csv")
mydata <- filter(mydata, tips != "Riebrev" & tips != "Leiolep")

# Read in the tree
mytree <- read.tree("chromEvol/tree_for_chromevol_exclude")

# Reorder data so it is the same order as the tree tip labels
mydata <- mydata[match(mytree$tip.label, mydata$tips), ]
# Replace tips with Binomials
mytree$tip.label <- mydata$Binomial

#------------------------------------------------
# Read in and process the ancestral states data
#------------------------------------------------
# Create anc states
ancs <- processAncStates("ChromoSSE/output/ChromoSSE_exclude_n18_final.tree", 
                     labels_as_numbers = TRUE)

# Change tip labels to species names
ancs@phylo$tip.label <- mytree$tip.label

# Get node age (node 92)
library(paleotree)
dateNodes(ancs@phylo)

#------------------------------------------------
# With tip labels
#------------------------------------------------

p2 <- 
  plotAncStatesMAP(ancs, 
                   tip_labels = TRUE,
                   #tip_states = TRUE,
                   #tip_states_size = 5,
                   #tip_labels_states = TRUE,
                   #tip_labels_states_size = 2,
                   #tip_labels_states_offset = 0,
                   node_labels_as = "state",
                   node_labels_centered = TRUE,
                   node_labels_offset = 0,
                   #node_size_as = "state_posterior",
                   alpha = 0.3)

p2 <- p2 +
  theme(legend.position = "none")

# let's add a time scale for key dates
root_age = 43.73925
p2 <- p2 + scale_x_continuous(breaks=c(root_age-33.9, root_age-30, root_age-23, root_age-18, root_age-14, root_age-10))
p2 <- p2 + theme(panel.grid.major = element_line(colour = "black", linewidth = 0.5, linetype = "dashed"),
                 panel.grid.major.y = element_blank()) 

p2

# Add biogeography
# Extract biogeography for the plotting
biog <-
  mydata %>%
  select(Realm) %>%
  as.data.frame()
# Add rownames as species names  
rownames(biog) <- mydata$Binomial

# Add to plot
p3 <- 
  gheatmap(p2, biog, offset = 12, width =.05, colnames = FALSE) +
  scale_fill_viridis_d(na.value = "grey90", option = "plasma", name = "")

# Look at this
#p3

ggsave(p3, file = "outputs/biogeography-figure.png", height = 10, width = 10, dpi = 600)

# Add letters in inkscape