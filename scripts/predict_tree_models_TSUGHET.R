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
library(climr)
library(sf)
library(raster)#masks dplyr select!!
library(ranger)

#load trained mods----
#set path to spp model outputs 
folder_path <- "outputs/ranger/RFregression_classes/TSUGHET"  
# List all RData files in the folder
Hwmods <- list.files(path = folder_path, pattern = "\\.RData$", full.names = TRUE)

# Load  & rename all mod files
for (i in seq_along(Hwmods)) {
  # Create a temporary environment to load the model
  temp_env <- new.env()
  
  # Load the model into the temporary environment
  load(Hwmods[i], envir = temp_env)  
  
  # Extract the new model name from the file path
  new_model_name <- stringr::str_match(Hwmods[i], "TSUGHET/(.*?)\\.RData")
  
  if (!is.na(new_model_name[,2])) {  # Check if the name extraction was successful
    # Get the name of mod in the temporary environment
    loaded_objects <- ls(envir = temp_env)
    
    model_name <- loaded_objects[1]  
    
    # Rename the model in the global environment
    assign(new_model_name[,2], get(model_name, envir = temp_env), envir = .GlobalEnv)
    
    # Optionally remove the old model if needed
    # rm(list = model_name, envir = .GlobalEnv)  # Uncomment if you want to clean up
  } else {
    message(paste("Could not extract name from:", Hwmods[i]))
  }
}

gc()

#predict----
#predict models over the full climate space of BC for each edatope 
#load downscaled climate data from PRISM DEM (generated in climr_getdata_projections.R) 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs (IN THE OCEAN) with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#add species - start with Hw
clim.bc$Species<-"TSUGHET"

#subset climate predictors 
clim<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim.bc<-dplyr::select(clim.bc, clim)  

#add loop through these.... 
#add edaphic info - start with Zonal  
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#make predictions
predsHwC4<-predict(object = RF_TSUGHET_C34,data =clim.bc, type='response')

#update edaphic info   
clim.bc$NutrientRegime_clean<-"A"   
clim.bc$MoistureRegime_clean<-"2" 

#make predictions
predsHwA2<-predict(object = RF_TSUGHET_AB12,data =clim.bc, type='response')

#update edaphic info   
clim.bc$NutrientRegime_clean<-"E"   
clim.bc$MoistureRegime_clean<-"6" 

#make predictions
predsHwE6<-predict(object = RF_TSUGHET_DE56,data =clim.bc, type='response')


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

predsHwC4.df<-process_predictions(predsHwC4)
predsHwA2.df<-process_predictions(predsHwA2)
predsHwE6.df<-process_predictions(predsHwE6)

##add lat, long back in
load(file="data/clim.bc2k.RData") #2km
predsHwA2.df<-left_join(predsHwA2.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsHwC4.df<-left_join(predsHwC4.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsHwE6.df<-left_join(predsHwE6.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)


#plot predictions---- 
library(tidyterra)

#load BC boundary for masking 
bcboundary<-bcmaps::bc_bound_hres()
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84 

  # Convert the data frame to a SpatRaster
  predsHwC4.r <- rast(predsHwC4.df, crs = "EPSG:4326")
  predsHwA2.r <- rast(predsHwA2.df, crs = "EPSG:4326")
  predsHwE6.r <- rast(predsHwE6.df, crs = "EPSG:4326")
  
  # Mask areas not in the BC boundary
  predsHwC4.r <- mask(predsHwC4.r, bcbound.reproj)
  predsHwA2.r <- mask(predsHwA2.r, bcbound.reproj)
  predsHwE6.r <- mask(predsHwE6.r, bcbound.reproj)
  
  # Convert the SpatRaster to a tidy data frame format for ggplot
  predsHwC4_tidy <- as_tibble(predsHwC4.r, xy = TRUE)
  predsHwA2_tidy <- as_tibble(predsHwA2.r, xy = TRUE)
  predsHwE6_tidy <- as_tibble(predsHwE6.r, xy = TRUE)
  
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
Hw<-subset(tree_dat , Species=='TSUGHET') 

# Plot using ggplot2
#plots<-
##  ggplot()+ geom_raster(data = predsHwC4_tidy, aes(x = x, y = y, fill = round(ypred)
#  geom_raster(data = bcbound.reproj, aes(x = x, y = y))
#  geom_point(data = Hw, aes(x = Longitude, y = Latitude, fill= TotalA_class), shape=21, size=3) +
#  scale_fill_grass_c(palette = 'water', direction = -1,
#                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%",">50%"))+
#  labs(title = "Hw plot data") + 
#  theme_classic() 

  
c4<-ggplot() +
    geom_raster(data = predsHwC4_tidy, aes(x = x, y = y, fill = round(ypred))) +
    scale_fill_grass_c(palette = 'water', direction = -1, #breaks=sort(unique(tree_dat$TotalA_class)),
                       name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%",">50%")) +
    geom_point(data = subset(Hw,NutrientRegime_clean=="C"& MoistureRegime_clean=="4"), 
               aes(x = Longitude, y = Latitude), color = "gold", size = 0.1) +
    labs(title = "Hw C4") + 
    theme_classic() #+ 
    #theme(legend.title = "abundance class", legend.text = )  # Turn off the legend if you want
a2<-ggplot() +
  geom_raster(data = predsHwA2_tidy, aes(x = x, y = y, fill = round(ypred))) +
  scale_fill_grass_c(palette = 'water', direction = -1, #breaks=sort(unique(tree_dat$TotalA_class)),
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%")) +
  geom_point(data = subset(Hw,NutrientRegime_clean=="A"& MoistureRegime_clean=="2"), 
             aes(x = Longitude, y = Latitude), color = "gold", size = 0.1) +
  labs(title = "Hw A2") + 
  theme_classic() 

e6<-ggplot() +
  geom_raster(data = predsHwE6_tidy, aes(x = x, y = y, fill = round(ypred))) +
  scale_fill_grass_c(palette = 'water', direction = -1, #breaks=sort(unique(tree_dat_sub$TotalA_class)),
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%")) +
  geom_point(data = subset(Hw,NutrientRegime_clean=="E"& MoistureRegime_clean=="6"), 
             aes(x = Longitude, y = Latitude), color = "gold", size = 0.1) +
  labs(title = "Hw E6") + 
  theme_classic() 


ggpubr::ggarrange(c4, a2, e6)

#extra code----
  
  #create all edatopic combinations for predictions
  edatopes<-expand.grid(NutrientRegime_clean= unique(tree_dat_wzeros$NutrientRegime_clean), 
                        MoistureRegime_clean= unique(tree_dat_wzeros$MoistureRegime_clean))%>%
    subset(MoistureRegime_clean!="1"& MoistureRegime_clean!="3"&MoistureRegime_clean!="5"&MoistureRegime_clean!="8"&
             NutrientRegime_clean!="B"& NutrientRegime_clean!="D"& NutrientRegime_clean!="B" & NutrientRegime_clean!="F")
  
  # Sample data frames
  
  # Create empty cols to fill with info 
  clim.bc$NutrientRegime_clean<-""
  clim.bc$MoistureRegime_clean<-""
  
  
  edatopes<-edatopes[c(1,5,9), ] #start with A2, C4, E6 only (like CCISS spatial)
  
  