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
#multi edatope model
load(file="outputs/ranger/RFregression_classes/RF_PSEUMEN.RData")
RF_PSEUMEN<-RF
rm(RF)

#predict----
#predict models over the full climate space of BC for each edatope 
#load downscaled climate data from PRISM DEM (generated in climr_getdata_projections.R) 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs (IN THE OCEAN) with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#add species
clim.bc$Species<-"PSEUMEN"

#subset climate predictors 
clim<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
        "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "Species")
clim.bc<-dplyr::select(clim.bc, clim)  


#create all edatopic combinations for predictions
load(file="data/tree_data_cleaned_wzeros.Rdata") 
tree_dat<-tree_dat_wzeros
rm(tree_dat_wzeros)

tree_dat<-mutate(tree_dat, TotalA_class=case_when(
  TotalA>50 & TotalA<=100~5,
  TotalA>25 & TotalA<=50~4,
  TotalA>=5 & TotalA<=25~3,
  TotalA>=1 & TotalA<5~2,
  TotalA>0 & TotalA<1~1,
  TotalA==0 ~0, 
  TRUE~NA))
Fd<-subset(tree_dat , Species=='PSEUMEN') 

edatopes<-expand.grid(NutrientRegime_clean= unique(Fd$NutrientRegime_clean), 
                      MoistureRegime_clean= unique(Fd$MoistureRegime_clean))

#edatopes<-edatopes[c(2,5,9), ] #start with B2, C4, E6 only (like CCISS spatial)

# Create empty cols to fill with info 
clim.bc$NutrientRegime_clean<-""
clim.bc$MoistureRegime_clean<-""

#loop to select each edatope and fit predictions model 








#add edaphic info - start with Zonal  
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#make predictions
predsFdC4<-predict(object = RF_PSEUMEN,data =clim.bc, type='response') 

#update edaphic info   
clim.bc$NutrientRegime_clean<-"B"   
clim.bc$MoistureRegime_clean<-"2" 

#make predictions
predsFdB2<-predict(object = RF_PSEUMEN,data =clim.bc, type='response') #multi edatope

#update edaphic info   
clim.bc$NutrientRegime_clean<-"E"   
clim.bc$MoistureRegime_clean<-"6" 

#make predictions
predsFdE6<-predict(object = RF_PSEUMEN,data =clim.bc, type='response') #multi edatope


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
predsFdB2.df<-process_predictions(predsFdB2)
predsFdE6.df<-process_predictions(predsFdE6)

##add lat, long back in
load(file="data/clim.bc2k.RData") #2km
predsFdB2.df<-left_join(predsFdB2.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsFdC4.df<-left_join(predsFdC4.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsFdE6.df<-left_join(predsFdE6.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)

#plot predictions---- 
library(tidyterra)

#load BC boundary for masking 
bcboundary<-bcmaps::bc_bound_hres()
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84 
#plot(st_geometry(bcbound.reproj))

# Convert the data frame to a SpatRaster
predsFdC4.r <- rast(predsFdC4.df, crs = "EPSG:4326")
predsFdB2.r <- rast(predsFdB2.df, crs = "EPSG:4326")
predsFdE6.r <- rast(predsFdE6.df, crs = "EPSG:4326")

# Mask areas not in the BC boundary
predsFdC4.r <- mask(predsFdC4.r, bcbound.reproj)
predsFdB2.r <- mask(predsFdB2.r, bcbound.reproj)
predsFdE6.r <- mask(predsFdE6.r, bcbound.reproj)

# Convert the SpatRaster to a tidy data frame format for ggplot
dummy<-data.frame(x=NA, y=NA, ypred=4)#add to preds if needed to have same number of classes for all edatopes 
dummy2<-data.frame(x=NA, y=NA, ypred=5)#add to preds if needed to have same number of classes for all edatopes 

predsFdC4_tidy <- as_tibble(predsFdC4.r, xy = TRUE)
predsFdB2_tidy <- as_tibble(predsFdB2.r, xy = TRUE)
predsFdE6_tidy <- as_tibble(predsFdE6.r, xy = TRUE)%>%rbind(.,dummy, dummy2)


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
B2<-ggplot() +
  geom_raster(data = predsFdB2_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Fd B2") +   xlab(" ") + ylab(" ")+ 
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
pdf(file = "outputs/ranger/RFregression_classes/pred_maps/PSEUMEN_preds_B2C4E6.pdf", width = 12, height = 8)
ggpubr::ggarrange(plots, c4,l1, B2, e6, l2)
dev.off()

#extract preds from plot locations for cross validation
Fdplots_sf = st_as_sf(Fd, coords = c("Longitude", "Latitude"), 
                      crs = 4326, agr = "constant")
CV_FdC4<-terra::extract(predsFdC4.r, Fdplots_sf)%>%cbind(., select(Fd, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="C" & MoistureRegime_clean=="4")

CV_FdB2<-terra::extract(predsFdB2.r, Fdplots_sf)%>%cbind(., select(Fd, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="B" & MoistureRegime_clean=="2")

CV_FdE6<-terra::extract(predsFdE6.r, Fdplots_sf)%>%cbind(., select(Fd, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="E" & MoistureRegime_clean=="6")

cor.test(CV_FdC4$ypred, CV_FdC4$TotalA_class) #80% 
ggplot(data=CV_FdC4, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Fd C4 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")

cor.test(CV_FdB2$ypred, CV_FdB2$TotalA_class) #80% 
ggplot(data=CV_FdB2, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Fd B2 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")

cor.test(CV_FdE6$ypred, CV_FdE6$TotalA_class) #78% #very few obs
ggplot(data=CV_FdE6, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Fd E6 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")

#single edatope models----
folder_path <- "outputs/ranger/RFregression_classes/Single_edatope/PSEUMEN"  
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

#make predictions
predsFdC4<-predict(object = RF_PSEUMEN_C34,data =clim.bc, type='response') #single edatope
predsFdB2<-predict(object = RF_PSEUMEN_AB12,data =clim.bc, type='response')#single edatope
predsFdE6<-predict(object = RF_PSEUMEN_DE56,data =clim.bc, type='response')#single edatope


