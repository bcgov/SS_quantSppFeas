#License info ----
#Copyright 2019 Province of British Columbia
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#libraries 
library(tidyverse)
library(terra)
library(data.table)
#library(climr)
library(sf)
library(raster)#masks dplyr select!!
library(ranger)


#load trained mods----
#set path to spp model outputs 
folder_path <- "outputs/ranger/RFregression_classes/PSEUMEN"  
# List all RData files in the folder
Fdmods <- list.files(path = folder_path, pattern = "\\.RData$", full.names = TRUE)

# Load  & rename all mod files
for (i in seq_along(Fdmods)) {
  # Create a temporary environment to load the model
  temp_env <- new.env()
  
  # Load the model into the temporary environment
  load(Fdmods[i], envir = temp_env)  
  
  # Extract the new model name from the file path
  new_model_name <- stringr::str_match(Fdmods[i], "PSEUMEN/(.*?)\\.RData")
  
  if (!is.na(new_model_name[,2])) {  # Check if the name extraction was successful
    # Get the name of mod in the temporary environment
    loaded_objects <- ls(envir = temp_env)
    
    model_name <- loaded_objects[1]  
    
    # Rename the model in the global environment
    assign(new_model_name[,2], get(model_name, envir = temp_env), envir = .GlobalEnv)
    
    # Optionally remove the old model if needed
    # rm(list = model_name, envir = .GlobalEnv)  # Uncomment if you want to clean up
  } else {
    message(paste("Could not extract name from:", Fdmods[i]))
  }
}

gc()

#predict----
#predict models over the full climate space of BC for each edatope 
#load downscaled climate data from PRISM DEM (generated in climr_getdata_projections.R) 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs (IN THE OCEAN) with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#add species - start with Fd
clim.bc$Species<-"PSEUMEN"

#subset climate predictors 
clim<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
        "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "Species")
clim.bc<-dplyr::select(clim.bc, clim)  

#add loop through these.... 
#add edaphic info - start with Zonal  
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#make predictions
predsFdC4<-predict(object = RF_PSEUMEN_C34,data =clim.bc, type='response')

#update edaphic info   
clim.bc$NutrientRegime_clean<-"A"   
clim.bc$MoistureRegime_clean<-"2" 

#make predictions
predsFdA2<-predict(object = RF_PSEUMEN_AB12,data =clim.bc, type='response')

#update edaphic info   
clim.bc$NutrientRegime_clean<-"E"   
clim.bc$MoistureRegime_clean<-"6" 

#make predictions
predsFdE6<-predict(object = RF_PSEUMEN_DE56,data =clim.bc, type='response')


#make into dfs
process_predictions <- function(preds) {
  # Convert predictions to a data frame
  preds_df <- as.data.frame(preds$predictions)
  
  # Create an ID column based on row names, converting to numeric
  preds_df$id <- as.numeric(as.character(row.names(preds_df)))
  
  # Rename columns to "ypred" and "id"
  colnames(preds_df) <- c("ypred", "id")
  
  return(preds_df)
}

predsFdC4.df<-process_predictions(predsFdC4)
predsFdA2.df<-process_predictions(predsFdA2)
predsFdE6.df<-process_predictions(predsFdE6)

##add lat, long back in
load(file="data/clim.bc2k.RData") #2km
predsFdA2.df<-left_join(predsFdA2.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsFdC4.df<-left_join(predsFdC4.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsFdE6.df<-left_join(predsFdE6.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)


#plot predictions---- 
library(tidyterra)

#load BC boundary for masking/ plotting  
bcboundary<-bcmaps::bc_bound_hres()
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84 
bcboundarylow<-bcmaps::bc_bound()
bcboundlow.reproj <- st_transform(bcboundarylow, st_crs(4326)) #reproject to wgs84 
#plot(st_geometry(bcbound.reproj))

# Convert the data frame to a SpatRaster
predsFdC4.r <- rast(predsFdC4.df, crs = "EPSG:4326")
predsFdA2.r <- rast(predsFdA2.df, crs = "EPSG:4326")
predsFdE6.r <- rast(predsFdE6.df, crs = "EPSG:4326")

# Mask areas not in the BC boundary
predsFdC4.r <- mask(predsFdC4.r, bcbound.reproj)
predsFdA2.r <- mask(predsFdA2.r, bcbound.reproj)
predsFdE6.r <- mask(predsFdE6.r, bcbound.reproj)


# Convert the SpatRaster to a tidy data frame format for ggplot
dummy<-data.frame(x=NA, y=NA, ypred=4)#add to preds if needed so all have 5 classes when plotting 
dummy2<-data.frame(x=NA, y=NA, ypred=5)#add to preds if needed so all have 5 classes when plotting 

predsFdC4_tidy <- as_tibble(predsFdC4.r, xy = TRUE)#%>%rbind(., dummy, dummy2)#already has all classes
predsFdA2_tidy <- as_tibble(predsFdA2.r, xy = TRUE)#%>%rbind(., dummy, dummy2)#already has all classes
predsFdE6_tidy <- as_tibble(predsFdE6.r, xy = TRUE)%>%rbind(., dummy, dummy2)



#load raw plot data for plotting 
load(file="data/tree_data_cleaned_wzeros.Rdata") 
tree_dat<-tree_dat_wzeros
tree_dat<-mutate(tree_dat, TotalA_class=case_when(
  TotalA>50 & TotalA<=100~5,
  TotalA>25 & TotalA<=50~4,
  TotalA>=5 & TotalA<=25~3,
  TotalA>=1 & TotalA<5~2,
  TotalA>0 & TotalA<1~1,
  TotalA==0 ~0, 
  TRUE~NA))
Fd<-subset(tree_dat , Species=='PSEUMEN') 

# Plot using ggplot2
#make dummy variable for first plot
predsFdC4_tidy$test<-0
predsFdC4_tidy<-mutate(predsFdC4_tidy, test=if_else(is.na(ypred),NA, 0))

plots<-ggplot()+ 
  geom_raster(data = predsFdC4_tidy, aes(x = x, y = y, fill = test)) +
  geom_point(data = Fd, aes(x = Longitude, y = Latitude, fill= TotalA_class), shape=21, size=3) +
  scale_fill_grass_c(palette = 'water', direction = -1,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%",">50%"))+
  labs(title = "Fd plot data")  + xlab(" ") + ylab(" ")+ 
   # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
c4<-ggplot() +
  geom_raster(data = predsFdC4_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Fd C4") +   xlab(" ") + ylab(" ")+ 
  #geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
l1<-ggpubr::get_legend(c4)
c4<-ggplot() +
  geom_raster(data = predsFdC4_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Fd C4") +   xlab(" ") + ylab(" ")+ 
  #  geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
a2<-ggplot() +
  geom_raster(data = predsFdA2_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
   labs(title = "Fd A2") +   xlab(" ") + ylab(" ")+ 
  #geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
   theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
e6<-ggplot() + 
  geom_raster(data = predsFdE6_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Fd E6") +   xlab(" ") + ylab(" ")+ 
 # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+
  theme_classic() + theme(axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank(), 
                          plot.margin = unit(c(2, 2, 2, 2),"cm"))
l2<-ggpubr::get_legend(e6)
e6<-ggplot() +
  geom_raster(data = predsFdE6_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
   scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Fd E6") +   xlab(" ") + ylab(" ")+ 
 # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none',axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
gc()

#combine into 4 panel 
pdf(file = "outputs/ranger/RFregression_classes/PSEUMEN/PSEUMEN_preds_A2C4E6.pdf", width = 12, height = 8)
ggpubr::ggarrange(plots, c4,l1, a2, e6, l2)
dev.off()
