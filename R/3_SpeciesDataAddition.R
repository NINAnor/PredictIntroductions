#---------------------------------------------------------------------#
# 3. Incorporate species specific data ####
#
# From now on we'll introduce species specific data. We need to get
# data for separate years though.
# 
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 3A. Narrow to Species-Specific Columms ####
#
#---------------------------------------------------------------------#

### We need to narrow down our data to the relevant species columns.
all_data <- readRDS("./Data/all_data.RDS")

# Create a vector to get rid of all the other species
species_to_remove <- species_list[species_list!=focal_species]

# Get only the names of the columns we ewant to keep
species_colnames <- colnames(all_data)[grepl(gsub(" ","_",focal_species),colnames(all_data))]
species_colnames_toKeep <- colnames(all_data)[!colnames(all_data) %in% 
                                                  grep(paste0(gsub(" ","_",species_to_remove), collapse = "|"), 
                                                       colnames(all_data), value = T)]

species_data <- all_data[,species_colnames_toKeep]

# Change names of species specific columns to a more friendly format
colnames(species_data)[colnames(species_data) %in% species_colnames] <- c("native","presence","introduced")

species_introSites <- species_data[species_data$introduced == 1,"no_vatn_lnr"]


#---------------------------------------------------------------------#
# 3B. Segment Data to Get Time-Dependent Variables ####
#
# We need to get variables which depend on the year (ie. what was
# the nearest population agt the time of introduction, not right now).
#
#---------------------------------------------------------------------#


# Need to read in original raw data to get years of introduction
raw_data <- readRDS("./Data/allSpecies_GBIFDownload_Edited.RDS")
raw_data$year <- ifelse(raw_data$datasetKey == "e306fa70-381e-4330-8e68-1f447b46a850", 1918, raw_data$year)

raw_data$scientificNameShort <- word(raw_data$scientificName, 1,2, sep=" ")

# Get earliest year with observation for every introduction.
get_years <- raw_data %>%
  as.data.frame() %>%
  filter(scientificNameShort == focal_species & vatnLnr %in% species_introSites & !is.na(year)) %>%
  dplyr::select(year, vatnLnr) %>%
  group_by(vatnLnr) %>%
  summarize(firstYear = min(year,na.rm=TRUE))


# We split the data into three.

# Firstly we produce a data frame of our introductions, excluding everything in the native range
# and everything with an absence, plus everything with no identifiable year.
species_intro_data <- species_data %>% 
  filter(introduced == 1 & no_vatn_lnr %in% get_years$vatnLnr)
species_intro_data$year <- ifelse(species_intro_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)

# We also produce a dataset with all presences, minus any introductions that have no identifiable
# year. These won't be in the model, but they'll be used to calculate nearest population
# and number of nearby populations.
species_prese_data <- species_data %>%
  filter(presence == 1 & 
           !(introduced ==1 & !(no_vatn_lnr %in% get_years$vatnLnr)))
species_prese_data$year <- ifelse(species_prese_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)

# We also need to produce a series of absences.
species_absen_data <- species_data %>%
  filter(native == 0 & presence == 0)
  
# We calculate nearby populations for the absences first. We take all current presences into account here,
# there's no need to put in a specific year.
# Calculate nearest population first.
nn_nearest_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=1)
species_absen_data$dist_n_pop <- nn_nearest_abs$nn.dist[,1]

# THen get all nearby populations.
nn_nearby_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=70)

# File this down to all populations within 5000m.
number_nearby_pop_5km <- nn_nearby_abs$nn.dist
number_nearby_pop_5km[number_nearby_pop_5km < 5000] <- 1
number_nearby_pop_5km[number_nearby_pop_5km >= 5000] <- 0
number_nearby_pop_5km_vec <- rowSums(number_nearby_pop_5km)

# And 10km, and 20km
number_nearby_pop_10km <- nn_nearby_abs$nn.dist
number_nearby_pop_10km[number_nearby_pop_10km < 10000] <- 1
number_nearby_pop_10km[number_nearby_pop_10km >= 10000] <- 0
number_nearby_pop_10km_vec <- rowSums(number_nearby_pop_10km)

number_nearby_pop_20km <- nn_nearby_abs$nn.dist
number_nearby_pop_20km[number_nearby_pop_20km < 20000] <- 1
number_nearby_pop_20km[number_nearby_pop_20km >= 20000] <- 0
number_nearby_pop_20km_vec <- rowSums(number_nearby_pop_20km)

# Put this into a data frame which we will then add more to later on.
spatial_data <- data.frame(species_absen_data$no_vatn_lnr, nn_nearest_abs$nn.dist, 
                           number_nearby_pop_5km_vec, number_nearby_pop_10km_vec, number_nearby_pop_20km_vec)
colnames(spatial_data) <- c("no_vatn_lnr","dist_n_pop","no_n_pop_5km", "no_n_pop_10km", "no_n_pop_20km")

# Now we cycle through year by year. All presences in native range are assumed to have been there ebfore the
# first time-step.
time_steps <- sort(unique(species_intro_data$year)) 
for (i in 1:length(time_steps)) {
  year_step <- time_steps[i]
  species_intro_historic <- species_intro_data %>% filter(year == year_step)
  species_presence_historic <- species_prese_data %>% filter(year < year_step)
  
  # Special case for first time step for Rainbow Trout, since they have no historic range in Norway.
  if (nrow(species_presence_historic) == 0) {
    spatial_data_timeStep <- data.frame(species_intro_historic$no_vatn_lnr, 9999999, 0, 0, 0)
    colnames(spatial_data_timeStep) <- c("no_vatn_lnr","dist_n_pop","no_n_pop_5km", "no_n_pop_10km", "no_n_pop_20km")
  } else {
  
  nn_nearest <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=2)
  # Same as above, but need to make sure we don't include a lake as its own enarest population
  nearest_pop_vec <- ifelse(nn_nearest$nn.dist[,1]==0,nn_nearest$nn.dist[,2],nn_nearest$nn.dist[,1])
  
  k_nearby <- ifelse(nrow(species_presence_historic) < 100, nrow(species_presence_historic)-1, 100)
  
  nn_nearby <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=k_nearby)
  
  number_nearby_pop_5km <- nn_nearby$nn.dist
  number_nearby_pop_5km[number_nearby_pop_5km < 5000 & number_nearby_pop_5km > 0] <- 1
  number_nearby_pop_5km[number_nearby_pop_5km >= 5000] <- 0
  number_nearby_pop_5km_vec <- rowSums(number_nearby_pop_5km)
  
  # And 10km, and 20km
  number_nearby_pop_10km <- nn_nearby$nn.dist
  number_nearby_pop_10km[number_nearby_pop_10km < 10000] <- 1
  number_nearby_pop_10km[number_nearby_pop_10km >= 10000] <- 0
  number_nearby_pop_10km_vec <- rowSums(number_nearby_pop_10km)
  
  number_nearby_pop_20km <- nn_nearby$nn.dist
  number_nearby_pop_20km[number_nearby_pop_20km < 20000] <- 1
  number_nearby_pop_20km[number_nearby_pop_20km >= 20000] <- 0
  number_nearby_pop_20km_vec <- rowSums(number_nearby_pop_20km)
  
  spatial_data_timeStep <- data.frame(species_intro_historic$no_vatn_lnr, nearest_pop_vec, 
                                      number_nearby_pop_5km_vec, number_nearby_pop_10km_vec, number_nearby_pop_20km_vec)
  colnames(spatial_data_timeStep) <- c("no_vatn_lnr","dist_n_pop","no_n_pop_5km", "no_n_pop_10km", "no_n_pop_20km")
  }
  
  spatial_data <- rbind(spatial_data,spatial_data_timeStep)

  }

#---------------------------------------------------------------------#
# 3C. Format Data for Use in Model ####
#
# We need to get variables which depend on the year (ie. what was
# the nearest population agt the time of introduction, not right now).
#
#---------------------------------------------------------------------#

# Now that we have these variables, need to produce a data frame free of incomplete cases.
# HOWEVER, we still need lakes outside the native range to find total number of rpesences
# nearby, so we keep these
species_model_data <- merge(species_data,spatial_data, all.x=TRUE, by="no_vatn_lnr")

species_model_data <- species_model_data[complete.cases(species_model_data) | species_model_data$native == 1,]

# Find duplicates
# all_data %>% group_by(no_vatn_lnr) %>%
#   tally() %>%
#   filter(n>1)
# 
# all_data %>% filter(no_vatn_lnr == 39447)
# raw_data %>% filter(vatnLnr == 39447)

duplicated_vatn_Lnrs <- unique(species_model_data$no_vatn_lnr[duplicated(species_model_data$no_vatn_lnr)])
if (length(duplicated_vatn_Lnrs) != 0) {
  warning(paste0("You have rows with duplicated Norwegian lake numbers. Lakes ",
               paste(duplicated_vatn_Lnrs,collapse=", ")," are duplicated. Use the function 
               display_duplicates to show them. Rows contained in vector duplicated_vatn_lnr."))
}

display_duplicates <- function(duplicated_vatn_Lnrs) {
  species_model_data %>% filter(no_vatn_lnr %in% duplicated_vatn_Lnrs)
}

# We're now dealing with species-specific data. All the data will be transferred to one file. 
# Create that file if it hasn't already been made.
if (dir.exists(paste0("./Data/Species_Data")) == FALSE) {dir.create(paste0("./Data/Species_Data"))}


if (dir.exists(paste0("./Data/Species_Data/",gsub(' ','_',focal_species))) == FALSE
) {dir.create(paste0("./Data/Species_Data/",gsub(' ','_',focal_species)))}

saveRDS(species_model_data,paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/species_model_data.RDS"))

print("Finished species data for use in model, can be found in species file in data folder under
      species_model_data.RDS")




