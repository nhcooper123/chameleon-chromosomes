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

mykey <- '2745adcf154b3539fae004c16f1bd9a6ea3a5138cf466a00e33b0beac666f8c2'
# If copying this code note that you need **your own API key**
# ***Do not use my API key***
#---------------------------------------
# Read in data
#---------------------------------------
# Read in the maps
maps <- st_read("raw-data/CHAMELEONS/")

# Read in the taxonomic data from Uetz 
taxonomy <- read.csv("raw-data/Chameleonidae-Uetz.csv")

# Read in climate data
bioclim <- getData("worldclim", var = "bio", res = 10)

#BIO1	Annual Mean Temperature
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
# Valid resolutions are 0.5, 2.5, 5, and 10 (minutes of a degree)

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
for(i in 1:19){
  mapsX <- 
    raster::extract(x = climate[[i]],
                    y = maps2,
                    fun = mean,
                    df = TRUE, 
                    small = TRUE,
                    na.rm = TRUE)
  
  # Export
  output[ , i+1] <- mapsX[ , 1]
}

# Wrangle metadata to get means per species
climate <- 
  output %>%
  group_by(binomial) %>%
  summarise_all(list(mean, max, min)) 

#---------------------------------------
# Extra IUCN data
#---------------------------------------
# Extract IUCN categories from IUCN
categories <- rl_comp_groups('chameleons', key = mykey)

# Tidy the output
iucn <- categories$result %>%
  select(scientific_name, category) %>%
  rename(binomial = scientific_name)


habitats <- as.data.frame(NULL)

for(i in 1:201) {
  
  rl_habitats(as.character(cham1$binomial[i]), 
                          key = mykey, 
                          parse = TRUE)
  
  habitats[i, 2] <- xx$result$code
  habitats[i, 2] <- xx$result$code
  habitats[i, 2] <- xx$result$code
  habitats[i, 2] <- xx$result$code
  habitats[i, 2] <- xx$result$code
}
habitats


#---------------------------------------
# Combine data
#---------------------------------------


areas
iucn
