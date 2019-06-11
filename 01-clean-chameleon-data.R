# Clean the raw data and extract data for analyses
# This looks complicated but if the data are 
# updated it's easy to rerun! 
# June 2019

#---------------------------------------
# Load libraries and functions
#---------------------------------------
library(tidyverse)

# ID approximations function
approx.needed <- function(data){
  
  if(
    # Is there an approximate value?
    sum(data$measure_type == "approx") > 0 &
    # And no other values
    sum(data$measure_type != "approx") == 0 
  ){
    out <- "Yes"
  }
  else{
    out <- "No"
  }
  return(out)
}


#---------------------------------------
# Read in data and wrangle
#---------------------------------------
ds <- read_csv("raw-data/chameleon data entry - Data.csv")

# First subset out the sexual dimorphism/other descriptions
# as it's easier to deal with separately
dimorphism <- filter(ds, variable == "dimorphism" | variable == "other")

# Clean the raw data
cleaned <- 
  ds %>%
  # Remove sexual dimorphism data as it's easier to deal with separately
  # And exclude "other" data
  filter(variable != "dimorphism" & variable != "other") %>%
  # Convert remaining values to numeric
  mutate(value = as.numeric(value)) %>%
  # Convert inches into mm where needed
  mutate(value, value = if_else(units == "inches", value * 25.4, value)) %>%
  # Unite variable and measure type columns
  unite(variable_type, variable, measure_type, 
        sep = "_", remove = FALSE) %>%
  # Select only necessary columns
  select(binomial, variable_type, value, variable, measure_type) %>%
  # Remove duplicates
  distinct()

#------------------------------------------------  
# Are the approximate values needed?
# If there are other records then we can omit
# the approximations. Otherwise we need to 
# keep them
#------------------------------------------------ 
# Create empty data frame
approx <- data.frame(array(dim = c(100, 3))) 
names(approx) <- c("binomial", "variable", "approximation")

# Start counter 
z <-1 

# For each species and variable check whether 
# the approximate values are needed
for (i in 1:length(unique(cleaned$binomial))){
  # Select one species
  sp <- filter(cleaned, binomial == unique(cleaned$binomial)[i])
  
  for(j in 1:length(unique(sp$variable))){
    # Select one variable
    sp_var <- filter(sp, variable == unique(sp$variable)[j])
    
    # Set counter
    k <- j - 1
    
    approx[z+k, 1] <-unique(sp_var$binomial)
    approx[z+k, 2] <-unique(sp_var$variable)
    approx[z+k, 3] <- approx.needed(sp_var)
  }
  # Make sure next species goes onto correct line of output
  z <- z + k + 1 
}

# Save a record of what we used approximations with
approx_used <-
  approx %>%
  filter(approximation == "Yes")

write_csv(approx_used, path = "raw-data/data-using-approximations.csv")

# Combine with cleaned data and remove approx values where not needed
cleaned_approx_fixed <-
  cleaned %>%
  full_join(approx, by = c("binomial", "variable")) %>%
  unite(remove, c(measure_type, approximation), remove = FALSE) %>%
  filter(remove != "approx_No") %>%
  select(-remove, -approximation)

#------------------------------------------------ 
# Extract summary variables
# We need to deal with elevation and body size 
# separately
#------------------------------------------------ 
# Elevation - max, min
elevation_min <-
  cleaned_approx_fixed %>%
  filter(variable == "elevation" & measure_type != "max") %>%
  group_by(binomial) %>%
  summarise(elevation_min = min(value)) 
  
elevation_max <-
  cleaned_approx_fixed %>%
  filter(variable == "elevation" & measure_type != "min") %>%
  group_by(binomial) %>%
  summarise(elevation_max = max(value)) 

# Size variables - max only
max_size <-
  cleaned_approx_fixed %>%
  filter(variable != "elevation" & measure_type != "min" & measure_type != "mean") %>%
  group_by(binomial, variable) %>%
  summarise(max = max(value)) 

# Fill in the gaps for total_length where 
# values for males or females are available
# but were not recorded as species totals
max_size <- 
  max_size %>% 
  spread(variable, max) %>%
  # Create new column of max of males and females
  # This creates -Inf if either are NA (this is where the warnings come from)
  mutate(max_mf = max(total_length_female, total_length_male, na.rm = TRUE)) %>%
  # Replace -Inf with NA
  mutate(max_mf = na_if(max_mf, "-Inf")) %>%
  # If total length is not available but male/female max is, then us this
  mutate(total_length = if_else(is.na(total_length) & !is.na(max_mf), max_mf, total_length)) %>%
  # Remove max_mf
  select(-max_mf)

final <- 
  max_size %>%
  full_join(elevation_max) %>%
  full_join(elevation_min)

#----------------------------------------------
# Combine with other data
#----------------------------------------------

climate <- read_csv("data/chameleon-iucn-data.csv")
karyotype <- read_csv("raw-data/karyotype.csv")

all <- 
  karyotype %>%
  full_join(climate) %>%
  full_join(final) 

write_csv(all, path = "data/all-june12-csv")

allk <- all %>%
  filter(karyotype == "y") %>%
  select(binomial, range, total_length:elevation_min) %>%
  distinct() %>%
  gather(measure, value, range:elevation_min)

allk %>% group_by(measure) %>% summarise(sum(!is.na(value)))

ggplot(allk, aes(x = (value), fill = measure)) +
  geom_histogram(alpha = 0.5, bins = 10) +
  facet_wrap(~ measure, scales = "free") +
  theme_bw() +
  theme(legend.position = "none")
