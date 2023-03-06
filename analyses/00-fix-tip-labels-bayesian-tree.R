# Data wrangling with new tree
# March 2022
##############################################
# Match Bayesian tree with chromosome dataset
##############################################
# Load packages
library(ape)

# Read in tree
tree <- read.tree("data/Bayesian_chameleon_constrained_time.tre")

# Read in new tip label names
tips <- read.csv("data/tip-names.csv")

# Drop tips that do not have data
tree <- drop.tip(tree, tips$x[which(is.na(tips$tip_name))])
# Check number of tips = 137
length(tree$tip.label)

# Now change names in the tree so they match the data
for(i in 1:length(tree$tip.label)){
  # Identify 
  tree$tip.label[i] <- tips$tip_name[tips$x == tree$tip.label[i]]
}

# Save tree
write.tree(tree, file = "data/Bayesian-tree.tre")

##############################################
# Extract tree for "species" level analyses
# Can also do this directly with script 1
##############################################

# Read in the species level dataset
mydata <- read.csv("data/chromosome-data-species.csv")

# Remove taxa not in this dataset
tree2 <- drop.tip(tree, setdiff(tree$tip.label, mydata$tips))

# Save tree
write.tree(tree2, file = "data/Bayesian-tree-species.tre")
