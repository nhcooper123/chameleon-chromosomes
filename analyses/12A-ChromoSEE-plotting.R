# Plotting ancestral states on tree

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
mydata <- read_csv("data/chromosome-data-species.csv")
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

#------------------------------------------------
# Plot the ancestral states tree 
#------------------------------------------------
# For some reason certain packages also being loaded
# (I suspect ChromEvol) break this code. Restart R if this happens.
p1 <- 
  plotAncStatesMAP(ancs, 
                       tip_labels = FALSE,
                       tip_states = TRUE,
                       tip_states_size = 5,
                       tip_labels_states = TRUE,
                       tip_labels_states_size = 2,
                       tip_labels_states_offset = 0,
                       node_labels_as = "state",
                       node_labels_centered = TRUE,
                       node_labels_offset = 0,
                       node_size_as = "state_posterior",
                       tree_layout = "fan", 
                       alpha = 0.3)

#p1 <- p1 +
 # theme(legend.position = "none")

# annotate the genera
p1 <- p1 + geom_cladelab(node = 94, label = "Brookesia", align = TRUE, offset = 5.0, offset.text = 1, fontface = 3)
p1 <- p1 + geom_cladelab(node = 105, label = "Rhampholeon", align = TRUE, offset = 5.0, offset.text = 2, fontface = 3)
p1 <- p1 + geom_cladelab(node = 136, label = "Furcifer", align = FALSE, offset = 5.0, offset.text = 20, fontface = 3)
p1 <- p1 + geom_cladelab(node = 156, label = "Calumma", align = TRUE, offset = 5.0, offset.text = 1, fontface = 3)
p1 <- p1 + geom_cladelab(node = 116, label = "Bradypodion", align = FALSE, offset = 5.0, offset.text = 30, fontface = 3)
p1 <- p1 + geom_cladelab(node = 110, label = "Chamaeleo", align = FALSE, offset = 5.0, offset.text = 10, fontface = 3)
p1 <- p1 + geom_cladelab(node = 125, label = "Kinyongia", align = TRUE, offset = 5.0, offset.text = 3, fontface = 3)
p1 <- p1 + geom_cladelab(node = 127, label = "Trioceros", align = FALSE, offset = 5.0, offset.text = 3, fontface = 3)

# Look at it
p1

ggsave(p1, file = "outputs/ChromoSSE_plot.png", width = 12, height = 12, dpi = 900)

# Then modify in inkscape
