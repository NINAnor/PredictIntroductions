#---------------------------------------------------------------------#
# 7. Visualise results of the forecasting
#---------------------------------------------------------------------#

# This script visualises where species are likely to be found after a 2.1 degree increase in
# temperature. 



forecasts <- readRDS(file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts.RDS"))
lake_likelihoods <- apply(forecasts, 1, mean)

all_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data.RDS"))
analytics <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_analytics.RDS"))

raw_data <- all_data$raw_data

Uncertainty <- analytics$intervals$width
InitProbs <- analytics$intervals$mean

all_data_likelihoods <- cbind(raw_data, lake_likelihoods, Uncertainty,InitProbs)
all_data_likelihoods$UncertaintyLevel <- as.factor(cut(all_data_likelihoods$Uncertainty, c(0, 0.01, 0.05, 0.15, 1.2),
                                                labels = c("negligible","very low", "low","moderate")))
all_data_likelihoods$PredictedIntro <- as.factor(cut(all_data_likelihoods$lake_likelihoods, c(-0.01, 0.05, 0.2, 0.5, 1.01),
                                                       labels = c("very low", "low","moderate","high")))


summary(all_data_likelihoods)

saveRDS(all_data_likelihoods, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/introduction_likelihoods.RDS"))
write.csv(all_data_likelihoods, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/introduction_likelihoods.csv"))


# Produce a map of Norway to use
Norway<-getData("GADM", country="NO", level=0)
Norway1<-getData("GADM", country="NO", level=1)
Norway1_sub<-Norway1[!(Norway1@data$NAME_1 %in% c('Troms', 'Finnmark', 'Nordland')),]
par(mar=c(1,1,1,1))
plot(Norway1_sub)
Norway_df <- fortify(Norway)
Norway1_sub_df <- fortify(Norway1_sub)

# These are just elements to clear the grid
theme_opts <- list(theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()))

# Let's first plot the initial appearances

initial_appearances <-  ggplot(data=all_data_likelihoods[all_data_likelihoods$presence == 1,]) + 
  geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
  geom_hex(aes(x = all_data_likelihoods[all_data_likelihoods$presence == 1,"decimalLongitude"], 
               y = all_data_likelihoods[all_data_likelihoods$presence == 1,"decimalLatitude"]),
           binwidth=c(0.3,0.18)) + theme_opts +
  scale_fill_gradient2(mid="yellow", high="red", #colors in the scale
                       #midpoint=40,    #same midpoint for plots (mean of the range)
                       breaks=seq(0,axis_limit,axis_limit/4), #breaks in the scale bar
                       limits=c(0,axis_limit)) +
  labs(subtitle=paste0("Current distribution of ",focal_species))


# Now let's plot predicted appearances, using likelihood of introduction 
# over our set threshold
predicted_appearances <- ggplot(data=all_data_likelihoods[all_data_likelihoods$lake_likelihoods >= prob_threshold,]) + 
  geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
  geom_hex(aes(x = all_data_likelihoods[all_data_likelihoods$lake_likelihoods >= prob_threshold,"decimalLongitude"], 
               y = all_data_likelihoods[all_data_likelihoods$lake_likelihoods >= prob_threshold,"decimalLatitude"]),
           binwidth=c(0.3,0.18)) + theme_opts +
  scale_fill_gradient2(mid="yellow", high="red", #colors in the scale
                       #midpoint=40,    #same midpoint for plots (mean of the range)
                       breaks=seq(0,axis_limit,axis_limit/4), #breaks in the scale bar
                       limits=c(0,axis_limit)) +
  labs(subtitle=paste0("Distribution of ",focal_species," in 50 Years"))


# Let's first plot the initial appearances

Internal_Uncertainty <- ggplot(data=all_data_likelihoods[all_data_likelihoods$introduced == 1,]) + 
  geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
  geom_point(aes(x = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLongitude"], 
                 y = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLatitude"],
                 color=UncertaintyLevel)) + 
  scale_color_manual(values=c("black", "brown", "orange", "yellow")) +
  theme_opts +
  labs(subtitle=paste0("Current distribution of ",focal_species))


# Let's first plot the initial appearances

Internal_Prediction <- ggplot(data=all_data_likelihoods[all_data_likelihoods$introduced == 1,]) + 
  geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
  geom_point(aes(x = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLongitude"], 
               y = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLatitude"],
               color=PredictedIntro)) + 
  scale_color_manual(values=c("white", "yellow", "orange", "red")) +
               theme_opts +
  labs(subtitle=paste0("Current distribution of ",focal_species))

maps <- list(initial_appearances = initial_appearances, predicted_appearances = predicted_appearances,
             Internal_Prediction = Internal_Prediction, Internal_Uncertainty = Internal_Uncertainty)


saveRDS(maps, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/maps.RDS"))



