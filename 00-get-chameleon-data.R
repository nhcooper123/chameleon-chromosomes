# Code to collate data from IUCN
# Natalie Cooper
# May 2019

#---------------------------------------
# Load libraries
#---------------------------------------
library(sf)
library(tidyverse)
library(rredlist)
library(raster)
library(sp)
library(rgeos)

#IUCN_REDLIST_KEY = '2745adcf154b3539fae004c16f1bd9a6ea3a5138cf466a00e33b0beac666f8c2')

# mykey <- '2745adcf154b3539fae004c16f1bd9a6ea3a5138cf466a00e33b0beac666f8c2'
# If copying this code note that you need **your own API key**
# ***Do not use my API key*** I will find you.
#---------------------------------------
# Read in data
#---------------------------------------
# Read in the maps
maps <- st_read("raw-data/CHAMELEONS/")

# Read in the taxonomic data from Uetz 
taxonomy <- read.csv("raw-data/Chameleonidae-Uetz.csv")

# Read in climate data
# Valid resolutions are 0.5, 2.5, 5, and 10 (minutes of a degree)
# Although sata is made at 30 sec this is probably overly detailed
bioclim <- getData("worldclim", var = "bio", res = 5)

#BIO1	Annual Mean Temperature * 10
#BIO2	Mean Diurnal Range (Mean of monthly (max temp – min temp))
#BIO3	Isothermality (BIO2/BIO7) (* 100)
#BIO4	Temperature Seasonality (standard deviation *100)
#BIO5	Max Temperature of Warmest Month
#BIO6	Min Temperature of Coldest Month
#BIO7	Temperature Annual Range (BIO5-BIO6)
#BIO8	Mean Temperature of Wettest Quarter
#BIO9	Mean Temperature of Driest Quarter
#BIO10	Mean Temperature of Warmest Quarter
#BIO11	Mean Temperature of Coldest Quarter
#BIO12	Annual Precipitation
#BIO13	Precipitation of Wettest Month
#BIO14	Precipitation of Driest Month
#BIO15	Precipitation Seasonality (Coefficient of Variation)
#BIO16	Precipitation of Wettest Quarter
#BIO17	Precipitation of Driest Quarter
#BIO18	Precipitation of Warmest Quarter
#BIO19	Precipitation of Coldest Quarter


#-------------------------------
# Extract geographic range areas
#-------------------------------
# Remove geometry info to leave metadata only
cham <- maps
st_geometry(cham) <- NULL

# Wrangle metadata to sum the areas per species
areas <- cham %>%
  group_by(binomial) %>%
  summarise(range = sum(SHAPE_Area)) 

#---------------------------------------
# Extract climate data
# For each of the 19 bioclim variables extract
# mean for each polygon and add to maps2
#---------------------------------------
# Create empty data frame
output <- data.frame(array(NA, dim = (c(nrow(cham), 20))))
names(output) <- c("binomial",
                   "bio1", "bio2", "bio3", "bio4", "bio5",
                   "bio6", "bio7", "bio8", "bio9", "bio10",
                   "bio11", "bio12", "bio13", "bio14", "bio15",
                   "bio16", "bio17", "bio18", "bio19")

output$binomial <- cham$binomial

# Convert to Spatial polygons to use raster
maps2 <- as(maps, 'Spatial')

# Run for all 19 bioclim variables
for(i in 2:19){
  mapsX <- 
    raster::extract(x = bioclim[[i]],
                    y = maps2,
                    fun = mean,
                    df = TRUE, 
                    # Do we want to use values where the polygons are small?
                    small = TRUE,
                    # Remove NAs when calculating means?
                    na.rm = TRUE)
  
  # Export
  output[ , i+1] <- mapsX[ , 2]
}

# Wrangle metadata to get means per species
climate <- 
  output %>%
  group_by(binomial) %>%
  summarise_all(mean) 

write_csv(climate, path = "climate-interim.csv")

#---------------------------------------
# Extra IUCN data
#---------------------------------------
# Extract IUCN categories from IUCN
categories <- rl_comp_groups('chameleons', key = mykey)

# Tidy the output
iucn <- categories$result %>%
  select(scientific_name, category) %>%
  rename(binomial = scientific_name)

# Extract habitats - this has to be done one species at a time
# This is made more difficuly by some species appearing in > 1 habitat
habitats <- data.frame(array(dim = c(10, 6)))
names(habitats) <- c("binomial", "habitat_code", "habitat", "suitability", "major_importantance")
z <- 1

for(i in 1:length(unique(cham$binomial))) {
  
  xx <- rl_habitats(as.character(unique(cham$binomial)[i]), 
                    key = mykey, 
                    parse = TRUE)
  
  # Output one line at a time to deal with species in > 1 habitat
  for (j in 1:length(xx$result$code)){
    
    # Create counter to get results in correct slots
    k <- j - 1
  
    habitats[z+k, 1] <- xx$name
    habitats[z+k, 2] <- xx$result$code[[j]]
    habitats[z+k, 3] <- xx$result$habitat[[j]]
    habitats[z+k, 4] <- xx$result$suitability[[j]]
    habitats[z+k, 5] <- xx$result$season[[j]]
    habitats[z+k, 6] <- xx$result$majorimportance[[j]]
    
  }

  # Make sure next species goes onto correct line of output
  z <- z + k + 1
}

#---------------------------------------
# Combine data
#---------------------------------------

taxonomy
areas
iucn
climate
habitats 
