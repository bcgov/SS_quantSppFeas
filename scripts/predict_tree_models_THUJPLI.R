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
load(file="outputs/ranger/RFregression_classes/RF_THUJPLI.RData")
RF_THUJPLI<-RF
rm(RF)

#predict----
#predict models over the full climate space of BC for each edatope 
#load downscaled climate data from PRISM DEM (generated in climr_getdata_projections.R) 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs (IN THE OCEAN) with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#add species
clim.bc$Species<-"THUJPLI"

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
Cw<-subset(tree_dat , Species=='THUJPLI') 

edatopes<-expand.grid(NutrientRegime_clean= unique(Cw$NutrientRegime_clean), 
                      MoistureRegime_clean= unique(Cw$MoistureRegime_clean))

#edatopes<-edatopes[c(2,5,9), ] #start with B2, C4, E6 only (like CCISS spatial)

# Create empty cols to fill with info 
clim.bc$NutrientRegime_clean<-""
clim.bc$MoistureRegime_clean<-""

#loop to select each edatope and fit predictions model 

#add edaphic info - start with Zonal  
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#make predictions
predsCwC4<-predict(object = RF_THUJPLI,data =clim.bc, type='response') 

#update edaphic info   
clim.bc$NutrientRegime_clean<-"B"   
clim.bc$MoistureRegime_clean<-"2" 

#make predictions
predsCwB2<-predict(object = RF_THUJPLI,data =clim.bc, type='response') #multi edatope

#update edaphic info   
clim.bc$NutrientRegime_clean<-"E"   
clim.bc$MoistureRegime_clean<-"6" 

#make predictions
predsCwE6<-predict(object = RF_THUJPLI,data =clim.bc, type='response') #multi edatope


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

predsCwC4.df<-process_predictions(predsCwC4)
predsCwB2.df<-process_predictions(predsCwB2)
predsCwE6.df<-process_predictions(predsCwE6)

##add lat, long back in
load(file="data/clim.bc2k.RData") #2km
predsCwB2.df<-left_join(predsCwB2.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsCwC4.df<-left_join(predsCwC4.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)
predsCwE6.df<-left_join(predsCwE6.df, dplyr::select(clim.bc, id, lat, lon))%>%dplyr::select(lon, lat, ypred)%>%rename(x=lon, y=lat)

#plot predictions---- 
library(tidyterra)

#load BC boundary for masking 
bcboundary<-bcmaps::bc_bound_hres()
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84 
#plot(st_geometry(bcbound.reproj))

# Convert the data frame to a SpatRaster
predsCwC4.r <- rast(predsCwC4.df, crs = "EPSG:4326")
predsCwB2.r <- rast(predsCwB2.df, crs = "EPSG:4326")
predsCwE6.r <- rast(predsCwE6.df, crs = "EPSG:4326")

# Mask areas not in the BC boundary
predsCwC4.r <- mask(predsCwC4.r, bcbound.reproj)
predsCwB2.r <- mask(predsCwB2.r, bcbound.reproj)
predsCwE6.r <- mask(predsCwE6.r, bcbound.reproj)

# Convert the SpatRaster to a tidy data frame format for ggplot
dummy<-data.frame(x=NA, y=NA, ypred=4)#add to preds if needed to have same number of classes for all edatopes 
dummy2<-data.frame(x=NA, y=NA, ypred=5)#add to preds if needed to have same number of classes for all edatopes 

predsCwC4_tidy <- as_tibble(predsCwC4.r, xy = TRUE)
predsCwB2_tidy <- as_tibble(predsCwB2.r, xy = TRUE)%>%rbind(., dummy2)
predsCwE6_tidy <- as_tibble(predsCwE6.r, xy = TRUE)%>%rbind(., dummy2)


# Plot using ggplot2
#make dummy variable for first plot
predsCwC4_tidy$test<-0
predsCwC4_tidy<-mutate(predsCwC4_tidy, test=if_else(is.na(ypred),NA, 0))

plots<-ggplot()+ 
  geom_raster(data = predsCwC4_tidy, aes(x = x, y = y, fill = test)) +
  geom_point(data = Cw, aes(x = Longitude, y = Latitude, fill= TotalA_class), shape=21, size=3) +
  scale_fill_grass_c(palette = 'water', direction = -1,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%",">50%"))+
  labs(title = "Cw plot data")  + xlab(" ") + ylab(" ")+ 
  # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
c4<-ggplot() +
  geom_raster(data = predsCwC4_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Cw C4") +   xlab(" ") + ylab(" ")+ 
  #geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
l1<-ggpubr::get_legend(c4)
c4<-ggplot() +
  geom_raster(data = predsCwC4_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Cw C4") +   xlab(" ") + ylab(" ")+ 
  #  geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
B2<-ggplot() +
  geom_raster(data = predsCwB2_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Cw B2") +   xlab(" ") + ylab(" ")+ 
  #geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
e6<-ggplot() + 
  geom_raster(data = predsCwE6_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Cw E6") +   xlab(" ") + ylab(" ")+ 
  # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+
  theme_classic() + theme(axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank(), 
                          plot.margin = unit(c(2, 2, 2, 2),"cm"))
l2<-ggpubr::get_legend(e6)
e6<-ggplot() +
  geom_raster(data = predsCwE6_tidy, aes(x = x, y = y, fill = factor(round(ypred), levels=0:5))) +
  scale_fill_grass_d(palette = 'water', direction = -1, breaks=0:5,
                     name = "Abundance class", labels = c("absent", "<1%", "1-5%", "5-25%", "25-50%", ">50%")) +
  labs(title = "Cw E6") +   xlab(" ") + ylab(" ")+ 
  # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none',axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
gc()

#combine into 4 panel 
pdf(file = "outputs/ranger/RFregression_classes/pred_maps/THUJPLI_preds_B2C4E6.pdf", width = 12, height = 8)
ggpubr::ggarrange(plots, c4,l1, B2, e6, l2)
dev.off()

#k-fold CV----
library(blockCV)
#get plots in sf 
Cwplots_sf = st_as_sf(Cw, coords = c("Longitude", "Latitude"), 
                      crs = 4326, agr = "constant")





#single edatope models----
folder_path <- "outputs/ranger/RFregression_classes/Single_edatope/THUJPLI"  
# List all RData files in the folder
Cwmods <- list.files(path = folder_path, pattern = "\\.RData$", full.names = TRUE)

# Load  & rename all mod files
for (i in seq_along(Cwmods)) {
  # Create a temporary environment to load the model
  temp_env <- new.env()
  
  # Load the model into the temporary environment
  load(Cwmods[i], envir = temp_env)  
  
  # Extract the new model name from the file path
  new_model_name <- stringr::str_match(Cwmods[i], "THUJPLI/(.*?)\\.RData")
  
  if (!is.na(new_model_name[,2])) {  # Check if the name extraction was successful
    # Get the name of mod in the temporary environment
    loaded_objects <- ls(envir = temp_env)
    
    model_name <- loaded_objects[1]  
    
    # Rename the model in the global environment
    assign(new_model_name[,2], get(model_name, envir = temp_env), envir = .GlobalEnv)
    
    # Optionally remove the old model if needed
    # rm(list = model_name, envir = .GlobalEnv)  # Uncomment if you want to clean up
  } else {
    message(paste("Could not extract name from:", Cwmods[i]))
  }
}

#make predictions
predsCwC4<-predict(object = RF_THUJPLI_C34,data =clim.bc, type='response') #single edatope
predsCwB2<-predict(object = RF_THUJPLI_AB12,data =clim.bc, type='response')#single edatope
predsCwE6<-predict(object = RF_THUJPLI_DE56,data =clim.bc, type='response')#single edatope




#extract preds from plot locations for cross validation
Cwplots_sf = st_as_sf(Cw, coords = c("Longitude", "Latitude"), 
                      crs = 4326, agr = "constant")
CV_CwC4<-terra::extract(predsCwC4.r, Cwplots_sf)%>%cbind(., select(Cw, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="C" & MoistureRegime_clean=="4")

CV_CwB2<-terra::extract(predsCwB2.r, Cwplots_sf)%>%cbind(., select(Cw, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="B" & MoistureRegime_clean=="2")

CV_CwE6<-terra::extract(predsCwE6.r, Cwplots_sf)%>%cbind(., select(Cw, TotalA_class, NutrientRegime_clean, MoistureRegime_clean))%>%
  subset(., NutrientRegime_clean=="E" & MoistureRegime_clean=="6")

cor.test(CV_CwC4$ypred, CV_CwC4$TotalA_class) #77% 
cor.test(log(CV_CwC4$ypred+1), log(CV_CwC4$TotalA_class+1)) #78% 

ggplot(data=CV_CwC4, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Cw C4 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")

cor.test(CV_CwB2$ypred, CV_CwB2$TotalA_class) #67% 
ggplot(data=CV_CwB2, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Cw B2 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")

cor.test(CV_CwE6$ypred, CV_CwE6$TotalA_class) #88% #very few obs
ggplot(data=CV_CwE6, aes(x=as.factor(TotalA_class), y=ypred))+
  geom_boxplot()+ geom_jitter(alpha=0.75, shape=21)+
  theme_classic()+ggtitle("Cw E6 predictions")+ xlab("plot abundance class") + ylab("pred abundance class")



