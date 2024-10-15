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

#predict over full climate surface with trained OF model(s)----
#load climate all BC (from PRISM DEM) 
#load(file="data/clim.bc.RData") #800m 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#create functions to streamline this!!! ####

#load trained models for comparison---- 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")#trained on 70-30 split
RFord5<-RFord
load(file="outputs/ordinalForest/RFordmodel6.Rdata")#trained on 70-30 split
RFord6<-RFord
load(file="outputs/ordinalForest/RFordmodel7.Rdata")#trained on 70-30 split
RFord7<-RFord
rm(RFord)
load(file="outputs/ordinalForest/RFordmodel8.Rdata")#mod 7 trained on all data 
RFord8<-RFord
rm(RFord)
load(file="outputs/ordinalForest/RFordmodel9.Rdata")#mod 6 trained on all data
RFord9<-RFord
rm(RFord)

#add edaphic info - start with Zonal for all 
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#add species - start with Hw
clim.bc$Species<-"TSUGHET"

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim5<-dplyr::select(clim.bc, clim5)  
clim6<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
"PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim6<-dplyr::select(clim.bc, clim6)  
clim7<-c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
        "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim7<-dplyr::select(clim.bc, clim7)  
clim8<-c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
         "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim8<-dplyr::select(clim.bc, clim8) #same params as clim7  
clim9<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
   "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim9<-dplyr::select(clim.bc, clim9) #same params as clim6  

#make predictions for all BC 
library(ordinalForest)
predsm5c4Hw <-predict(object = RFord5, newdata = clim5)
save(predsm5c4Hw, file="outputs/ordinalForest/preds/predsm5c4Hw.RData")
predsm6c4Hw <-predict(object = RFord6, newdata = clim6)
save(predsm6c4Hw, file="outputs/ordinalForest/preds/predsm6c4Hw.RData")
predsm7c4Hw <-predict(object = RFord7, newdata = clim7)
save(predsm7c4Hw, file="outputs/ordinalForest/preds/predsm7c4Hw.RData")
predsm8c4Hw <-predict(object = RFord8, newdata = clim8)
save(predsm8c4Hw, file="outputs/ordinalForest/preds/predsm8c4Hw.RData")
save(predsm8c4Hw, file="outputs/ordinalForest/preds/predsc4Hw.RData")
predsm9c4Hw <-predict(object = RFord9, newdata = clim9)
save(predsm9c4Hw, file="outputs/ordinalForest/preds/predsm9c4Hw.RData")

#re-do for Cw
clim.bc$Species<-"THUJPLI"
clim8<-c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
         "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim8<-dplyr::select(clim.bc, clim8) #same params as clim7  
clim9<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim9<-dplyr::select(clim.bc, clim9) #same params as clim6  

predsm8c4Cw <-predict(object = RFord8, newdata = clim8)
save(predsm8c4Cw, file="outputs/ordinalForest/preds/predsm8c4Cw.RData")
predsm9c4Cw <-predict(object = RFord9, newdata = clim9)
save(predsm9c4Cw, file="outputs/ordinalForest/preds/predsm9c4Cw.RData")

#re-do for Fd
clim.bc$Species<-"PSEUMEN"
clim8<-c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
         "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim8<-dplyr::select(clim.bc, clim8) #same params as clim7  
clim9<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim9<-dplyr::select(clim.bc, clim9) #same params as clim6 
predsm8c4Fd <-predict(object = RFord8, newdata = clim8)
save(predsm8c4Fd, file="outputs/ordinalForest/preds/predsm8c4Fd.RData")
predsm9c4Fd <-predict(object = RFord9, newdata = clim9)
save(predsm9c4Fd, file="outputs/ordinalForest/preds/predsm9c4Fd.RData")



###Misc code----
#BEC plots with rank cover 
Hw1<-subset(Hw, cover_rank=="1")
points(x = Hw$Longitude, Hw$Latitude,
       col=c("blue", "yellow", "green")[Hw$cover_rank],
       pch=15,
       cex = 0.8)
legend("topright",   
       legend = levels(predsm5c4Hw$ypred), 
       # the legend
       fill =  c("forestgreen","lightgreen", "yellow", "green"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))
       
#mod 6
#load trained model 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-dplyr::select(clim.bc, clim5)  
clim5<-na.omit(clim5)  #remove NAs

#add edaphic info 
clim5$NutrientRegime_clean<-"C" #zonal  
clim5$MoistureRegime_clean<-"4" #zonal  
clim5c4<-clim5
clim5$NutrientRegime_clean<-"B" #dry poor 
clim5$MoistureRegime_clean<-"2" #dry poor
clim5b2<-clim5
clim5$NutrientRegime_clean<-"E" #wet rich 
clim5$MoistureRegime_clean<-"6" #wet rich 
clim5e6<-clim5

#add species 
clim5c4$Species<-"TSUGHET"
clim5b2$Species<-"TSUGHET"
clim5e6$Species<-"TSUGHET"

#make predictions for all BC 
rm(clim.bc)
gc()
predsm5c4Hw <- predict(object = RFord,newdata = clim5c4)
predsm5b2Hw <- predict(object = RFord,newdata = clim5b2)
predsm5e6Hw <- predict(object = RFord,newdata = clim5e6)



## populate the raster grid with prediction values
preds_map <- rast(dem.bc) # use the DEM as a template raster
preds_map[predsm5c4$id]] <- ds_out[PERIOD == "2001_2020", MAT] # populate the raster cells with the 2001-2020 mean annual temperature (MAT) values, using the `id` field as the link. 

plot(my_clim, main = "2001-2020 mean annual temperature (MAT)")


#mod 7
#load trained model 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-dplyr::select(clim.bc, clim5)  
clim5<-na.omit(clim5)  #remove NAs

#add edaphic info 
clim5$NutrientRegime_clean<-"C" #zonal  
clim5$MoistureRegime_clean<-"4" #zonal  
clim5c4<-clim5
clim5$NutrientRegime_clean<-"B" #dry poor 
clim5$MoistureRegime_clean<-"2" #dry poor
clim5b2<-clim5
clim5$NutrientRegime_clean<-"E" #wet rich 
clim5$MoistureRegime_clean<-"6" #wet rich 
clim5e6<-clim5

#add species 
clim5c4$Species<-"TSUGHET"
clim5b2$Species<-"TSUGHET"
clim5e6$Species<-"TSUGHET"

#make predictions for all BC 
rm(clim.bc)
gc()
predsm5c4Hw <- predict(object = RFord,newdata = clim5c4)
predsm5b2Hw <- predict(object = RFord,newdata = clim5c4)
predsm5e6Hw <- predict(object = RFord,newdata = clim5c4)

## populate the raster grid with prediction values
preds_map <- rast(dem.bc) # use the DEM as a template raster
preds_map[predsm5c4$id]] <- ds_out[PERIOD == "2001_2020", MAT] # populate the raster cells with the 2001-2020 mean annual temperature (MAT) values, using the `id` field as the link. 

plot(my_clim, main = "2001-2020 mean annual temperature (MAT)")












clim.bc.stack<-raster::stack(clim.bcr) #stack 

#save(clim.bc.stack, file="data/clim.bc.stack.RData")
plot(clim.bcr) #looks good 

#add nutrient & moisture regimes 
a<-raster(ncol=3241, nrow=1680)
a<-setValues(a, "A")
clim.bc.stack<-addLayer(clim.bc.stack, a)
b<-rast(ncol=3241, nrow=1680)
c<-rast(ncol=3241, nrow=1680)
d<-rast(ncol=3241, nrow=1680)
e<-rast(ncol=3241, nrow=1680)
f<-rast(ncol=3241, nrow=1680)

#predict onto full climate surface using trained OF models 
#mod5
load(file="data/spatial/clim.bc.stack.RData")



load(file="outputs/ordinalForest/RFordmodel5.Rdata")
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-select(clim.bc, clim5)  
clim5$NutrientRegime_clean<-"C" #use zonal here 
clim5$MoistureRegime_clean<-"3" #use zonal here 


spat_pred5<-raster::predict(clim5r, RFord, type="response")

rm(RFord)





---------------------
#add gps coords back in 
clim.bc<-left_join(clim.bc, my_grid)
clim.bc<-na.omit(clim.bc)  #remove NAs

clim.bc<-raster::rasterFromXYZ(clim.bc)

#lower resolution of latitute and longitude to 2 sig digits ~1km 
clim.bc$latlr<-round(clim.bc$lat, 3)
clim.bc$lonlr<-round(clim.bc$lon, 3)
#clim.bc<-group_by(clim.bc, latlr, lonlr)%>% mutate(elevlr=mean(elev)) #take mean elevs by GPS 

clim.bcx<-group_by(clim.bc, latlr, lonlr, elev)%>%   
  summarize(DD5_sp= mean(DD5_sp), Tmax_sm= mean(Tmax_sm), TD= mean(TD), PPT_sm=mean(PPT_sm))

#DDsub0_at= mean(DDsub0_at)
, "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
"PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFCw", "PAS", 
"CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp"))
gc()










#Predict for ALL species 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0
#select climate params from best model 
clim9<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFCw", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim9<-dplyr::select(clim.bc, clim9) 
#create dataframe with all spp
climspp<-expand.grid(clim9, Species=unique(tree_dat$Species))
#loop model over each spp dataset 
sppcodes <- sort(unique(ckun$geo))
for (k in codes) {
  name<-paste("data", k, sep="_")
  assign(name, subset(Data, (Data$geo==k)))
}

