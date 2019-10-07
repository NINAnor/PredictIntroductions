#---------------------------------------------------------------------#
# 6. Forecasting
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 6A. Get Scaling Parameters
#
# Variables within the model will change, but they have to be scaled 
# using the same mands and SDs as our original data, so we need to define
# those first.
#
#---------------------------------------------------------------------#


# Import data and draws and define components

if (use_buffered_model == TRUE) {
  model_output <- readRDS(paste0("./Data/",gsub(' ','_',focal_species),"/buffered_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/",gsub(' ','_',focal_species),"/buffered_model_data.RDS"))
} else {
  model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data.RDS"))
}

raw_data <- model_data_extra$raw_data

draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha

# Whatever data/model we are using for out parameters, we are now applying the parameters to all the data.
# That means we need to scale all the data using whatever means and SDs with which we scaled the data that
# was used for the moel.
model_data <- model_data_extra$raw_data %>%
  filter(native == 0)
#env_var <- c("HFP", "nearby_pops", "pop_dist", "area_km2", "perimeter_m", "distance_to_road", "eurolst_bio10")
data_4scaling <- model_data %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(dist_n_pop+1),
            SCI = log(((perimeter_m/1000)/(2*sqrt(pi*area_km2)))+1),
            HFP = HFP,
            n_pop = log(no_n_pop+1))

# Get the means and sds to work with
means_relData <- apply(data_4scaling, 2, mean)
SDs_relData <- apply(data_4scaling, 2, sd)

# Now we get our already scaled data
env_data <- model_data_extra$env_data

print("Data is scaled.")

#---------------------------------------------------------------------#
# 6B. Set Figures Beforehand
#
# We need to set up a bunch of parameters to use in the forecasting
# script. 
#
#---------------------------------------------------------------------#

### Get betas and alphas
calc_beta <- apply(as.matrix(calculate(beta,draws)),2,mean)
calc_alpha <- apply(as.matrix(calculate(alpha,draws)),2,mean)

# Need vector of where we have presences outside of the native range at the moment
presences <- model_data_extra$intro_data$introduced

# Get full index of nearby lakes that have a chance of colonisation, then narrow them down to lakes within 5000m
attempt <- as.matrix(env_data) %*% as.matrix(calc_beta)
eta <- sweep(attempt, 2, calc_alpha, "+")
expeta <- exp(eta)
init_probabilities <- expeta/(1+expeta)

nn_all <- get.knnx(raw_data[,c("utm_x","utm_y")],model_data[c("utm_x","utm_y")],k=100)
nn_all$nn.index[nn_all$nn.dist > 5000 | nn_all$nn.dist == 0] <- 0


periods <- list()

# These variables define actions in the loop
number_reps <- 5

# This defines in incremental increase over 50 years based on a 1 degree increase in temperature over 50 years. 10 represents a 1 degree increase.
temp_step <- 21/SDs_relData["temp"]

# These variables will be used to scale projections of numnber of populations nearby and distance to nearest population
mean_n_pop <- means_relData["n_pop"]
sd_n_pop <- SDs_relData["n_pop"]

sd_dist <- SDs_relData["pop_distL"]
mean_dist <- means_relData["pop_distL"]

introductions <- matrix(NA,nrow=nrow(env_data),ncol=1)
periods <- matrix(NA,nrow=nrow(env_data),ncol=n_loops)
# introduction_probs <- list()
# populations <- list()

print("Other stuff in place, starting forecast run.")


#---------------------------------------------------------------------#
# 6C. Run Forecasting Loop
#
#
#---------------------------------------------------------------------#


time1 <- Sys.time()
# data <- list()
for (s in 1:n_loops) {
  for (j in 1:number_reps) {
    # So, first step is to create predictions based purely off a rise in temperature. Because this is slightly different to what we'll do in the other steps, there needs to be an if clause.
    
    if (j==1) {
      # We create a table containing our new data
      
      tempIncrease <- rnorm(nrow(env_data), temp_step/5, 0.05)
      newData <- env_data
      newData[,"temp"] <- newData[,"temp"] + tempIncrease
      
      # Unfortunately because we're predicting for 600k + lakes, it's impossible to run the preditions all at one. So the loop below runs them in chunks of 10k lakes at a time.
      
      attempt <- as.matrix(newData) %*% as.matrix(calc_beta)
      eta <- greta::sweep(attempt, 2, calc_alpha, "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # # # Now figure out which first upstream lakes pike have spread to using the gradient you defined before. I wrote a function for this which has been loaded above.
      # introductions_by_dispersal <- introduction_by_dispersal(full_data[full_data$Esox_lucius==1,]$waterBodyID,200, connect_db)
      # new_presences <- ifelse(full_data$Esox_lucius ==1, 1, ifelse(threshold_var > threshold_int, 1, ifelse(full_data$waterBodyID %in% introductions_by_dispersal,1, 0)))
      # 
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      all_presences <- ifelse(presences ==1, 1, ifelse(new_presences == 1, 1, 0))
      
      newData2 <- newData
    } else {
      
      ### And now we move on to the second part of the loop, whereby we redefine those variables pertaining to closest population and numbers of close populations at each step
      
      # Introduce the new temperature thing first, it's the simplest
      tempIncrease <- rnorm(nrow(newData2), temp_step/5, 0.05)
      newData2[,"temp"] <- newData2[,"temp"] + tempIncrease
      
      # Now figure out the new distance to closest population
      pop_proximity <- cbind(all_presences, model_data[,c("utm_x","utm_y","decimalLatitude","decimalLongitude","locationID")])
      
      # Need to bind our new presences together with all the presences in the native range.
      data_presences_nonnative <- as.data.frame(pop_proximity %>% filter(all_presences == 1) %>% 
                                                  dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude))
      data_presences_native <- raw_data %>% filter(native == 1 & presence == 1) %>% 
                                                dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude)
      data_presences <- rbind(data_presences_nonnative, data_presences_native)
      
      data_all <- pop_proximity %>%
        dplyr::select(utm_x,utm_y,locationID, decimalLongitude, decimalLatitude) %>%
        distinct() %>%
        as.data.frame()
      
      # The get.knnx function returns the distance from a lake in table B to the k closest lakes in table A
      
      nn <- get.knnx(data_presences[c("utm_x","utm_y")],data_all[c("utm_x","utm_y")],k=2)
      
      # We then figure out the distance to the closest lake by taking the first column
      dist_to_closest_pop <- ifelse(nn$nn.dist[,1]==0,nn$nn.dist[,2],nn$nn.dist[,1])
      locationID <- data_all$locationID
      distance_data <- as.data.frame(cbind(dist_to_closest_pop,locationID))
      
      ordered_distances <- merge(model_data["locationID"],distance_data,all.x=TRUE,by="locationID")
      
      new_distances <- log(as.numeric(as.character(ordered_distances$dist_to_closest_pop))+1)
      
      # So now that we have the new measurements for closest population, these need to be scaled against those that we had for the initial population. So we subtract the mean and divide by the standard deviation.
      
      newData2[,"dist_n_popL"] <- (new_distances-mean_dist)/sd_dist
      
      # Now we need to get the new measurements for number of close populations
      
      all_presence_table <- raw_data["locationID"]
      all_presence_table$present <- ifelse(all_presence_table$locationID %in% data_presences$locationID, 1, 0)
      
      
      rowNumbers_withPresence <- which(all_presence_table$present==1)
      nn_inds <- nn_all$nn.index
      nn_inds[!(nn_inds %in% rowNumbers_withPresence)] <- 0
      nn_inds[nn_inds %in% rowNumbers_withPresence] <- 1
      number_nearby_pop_vec <- log(rowSums(nn_inds)+1)
      # 
      # Scale it like we did for the last variable
      
      newData2[,"no_n_popL"] <- (number_nearby_pop_vec-mean_n_pop)/sd_n_pop
      
      attempt <- as.matrix(newData2) %*% as.matrix(calc_beta)
      eta <- sweep(attempt, 2, calc_alpha, "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      
      # # Introduce upstream lakes again
      # introductions_by_dispersal2 <- introduction_by_dispersal(data_presences$waterBodyID,200, connect_db)
      # new_presences <- ifelse(pop_proximity$new_presences ==1, 1, ifelse(threshold_var > threshold_int, 1, ifelse(pop_proximity$waterBodyID %in% introductions_by_dispersal2,1, 0)))
      
      
      # Establish which lakes now have presences.
      all_presences <- ifelse(pop_proximity$all_presences ==1, 1, ifelse(new_presences == 1, 1, 0))
    }
  }
  time2 <- Sys.time()
  difftime_1 <- round(as.numeric(difftime(time2, time1,
                                          units = "mins")),3)
  if (s %% 5 == 0) {print(paste0("Run ", s, " of ",  n_loops, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*n_loops/s-difftime_1,5), " minutes left."))}
  periods[,s] <- all_presences
}


print("Finished forecasting.")


forecasts <- periods
saveRDS(forecasts, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts.RDS"))
summary(forecasts)
summary(rowSums(forecasts)/100)
