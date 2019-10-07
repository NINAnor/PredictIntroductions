### Overarchign script, tying everythign together

### Let's get all the libraries in

library(dplyr)
library(pool)
library(postGIStools)
library(sf)
library(getPass)
library(rgbif)
library(dplyr) # for data-wrangling
library(mapedit)
library(rio)
library(raster)
library(FNN)
library(readr)
library(ggplot2)
library(greta)
library(stringr)
source("ExtraFunctions.R")


### Script 1 downloads everything for you for all species. 

# Define species of interest

species_list <- c("Coregonus lavaretus", # define species of interest
                  "Esox lucius", 
                  "Rutilus rutilus",
                  "Scardinius erythrophthalmus", 
                  "Perca fluviatilis",
                  "Oncorhynchus mykiss")

# All the data you create will be transferred into a folder. Let's create that folder now.
if (dir.exists(paste0("./Data")) == FALSE
    ) {dir.create(paste0("./Data"))}

initiate_download <- FALSE            # Have you already initiated a download? If so, set to false.
GBIF_download <- FALSE                 # Have you already downloaded occurrence data from GBIF? If so, set to false.
download_lakes <- FALSE               # Have you already downloaded the lakes? If so, set to false.
get_all_lakes <-  FALSE               # Have you downloaded all lake data from NOFA? If so, set to false.
delete_north <-  TRUE                 # Do you want to use lakes north of Tronderlag? If so, set to true.
download_native_range <- FALSE        # Have you already downloaded species' native ranges? If so, set to false.

#   dist_threshold sets the distance between a point and its designated lake which
#   is acceptable.
dist_threshold <- 50

source("./R/1_GetIntroductions.R")

# Next step is a fairly short script to add in all environmental covariates which can be 
# calculated regardless of species.


size_threshold <- 0.02

HFP_download <- TRUE                # Set this to false if you already have raster data downloaded in your data folder.

source("./R/2_BioticDataAddition.R")

# From now on everything becomes species specific.

focal_species <- species_list[2]

# Now we run the preliminary model. Only thing we need to choose is what our size limit 
# on lakes will be.

# Don't worry about the duplication if the only lake that has been duplicated is 39447.
# If there are others, let me know.

source("./R/3_SpeciesDataAddition.R")

# Define radius you want to measure nearby populations by
population_threshold <- 5

# This script will take a while, as it's running a model on up to 250,000 lakes. Grab a 
# coffee. Teach it to do algebra.
source("./R/4_FullModelConstruct.R")


# Next one gives you uncertainty from each lake, convergence diagnostic,
# beta intervals, and model deviance.
# 
use_buffered_model <- FALSE
source("./R/5_ModelAnalysis.R")

# Now we can run our simulations. need to define which model we want to use.
# n_loops simply tells us how many iterations we want to run. Be careful,
# because the run time can be enormous.
n_loops <- 100

source("./R/6_Looping.R")

# Last step is simply doing some data visualisation.

prob_threshold <- 0.5
axis_limit <- 150

source("./R/7_Visualisations.R")



maps$initial_appearances
maps$predicted_appearances


