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

#Load model predictions 
load("outputs/ordinalForest/preds/predsm5c4Hw.RData")
load("outputs/ordinalForest/preds/predsm6c4Hw.RData")
load("outputs/ordinalForest/preds/predsm7c4Hw.RData")
load("outputs/ordinalForest/preds/predsm8c4Hw.RData")
load("outputs/ordinalForest/preds/predsm8c4Cw.RData")
load("outputs/ordinalForest/preds/predsm8c4Fd.RData")
load("outputs/ordinalForest/preds/predsm9c4Hw.RData")
load("outputs/ordinalForest/preds/predsm9c4Cw.RData")
load("outputs/ordinalForest/preds/predsm9c4Fd.RData")

#create functions to streamline this!!! ####

#make into dfs
predsm5c4Hw<-as.data.frame(predsm5c4Hw$ypred) 
predsm5c4Hw$id<-as.numeric(as.character(row.names(predsm5c4Hw)))
colnames(predsm5c4Hw)<- c("ypred", "id")

predsm6c4Hw<-as.data.frame(predsm6c4Hw$ypred) 
predsm6c4Hw$id<-as.numeric(as.character(row.names(predsm6c4Hw)))
colnames(predsm6c4Hw)<- c("ypred", "id")

predsm7c4Hw<-as.data.frame(predsm7c4Hw$ypred) 
predsm7c4Hw$id<-as.numeric(as.character(row.names(predsm7c4Hw)))
colnames(predsm7c4Hw)<- c("ypred", "id")

predsm8c4Hw<-as.data.frame(predsm8c4Hw$ypred) 
predsm8c4Hw$id<-as.numeric(as.character(row.names(predsm8c4Hw)))
colnames(predsm8c4Hw)<- c("ypred", "id")

predsm9c4Hw<-as.data.frame(predsm9c4Hw$ypred) 
predsm9c4Hw$id<-as.numeric(as.character(row.names(predsm9c4Hw)))
colnames(predsm9c4Hw)<- c("ypred", "id")

predsm8c4Cw<-as.data.frame(predsm8c4Cw$ypred) 
predsm8c4Cw$id<-as.numeric(as.character(row.names(predsm8c4Cw)))
colnames(predsm8c4Cw)<- c("ypred", "id")

predsm9c4Cw<-as.data.frame(predsm9c4Cw$ypred) 
predsm9c4Cw$id<-as.numeric(as.character(row.names(predsm9c4Cw)))
colnames(predsm9c4Cw)<- c("ypred", "id")

predsm8c4Fd<-as.data.frame(predsm8c4Fd$ypred) 
predsm8c4Fd$id<-as.numeric(as.character(row.names(predsm8c4Fd)))
colnames(predsm8c4Fd)<- c("ypred", "id")

predsm9c4Fd<-as.data.frame(predsm9c4Fd$ypred) 
predsm9c4Fd$id<-as.numeric(as.character(row.names(predsm9c4Fd)))
colnames(predsm9c4Fd)<- c("ypred", "id")

##add lat, long back in
load(file="data/clim.bc2k.RData") #2km
predsm5c4Hw<-left_join(predsm5c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm6c4Hw<-left_join(predsm6c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm7c4Hw<-left_join(predsm7c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm8c4Hw<-left_join(predsm8c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm9c4Hw<-left_join(predsm9c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm8c4Cw<-left_join(predsm8c4Cw, dplyr::select(clim.bc, id, lat, lon))
predsm9c4Cw<-left_join(predsm9c4Cw, dplyr::select(clim.bc, id, lat, lon))
predsm8c4Fd<-left_join(predsm8c4Fd, dplyr::select(clim.bc, id, lat, lon))
predsm9c4Fd<-left_join(predsm9c4Fd, dplyr::select(clim.bc, id, lat, lon))

#rename for raster creation 
predsm5c4Hw<-dplyr::select(predsm5c4Hw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm6c4Hw<-dplyr::select(predsm6c4Hw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm7c4Hw<-dplyr::select(predsm7c4Hw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm8c4Hw<-dplyr::select(predsm8c4Hw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm9c4Hw<-dplyr::select(predsm9c4Hw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm8c4Cw<-dplyr::select(predsm8c4Cw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm9c4Cw<-dplyr::select(predsm9c4Cw, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm8c4Fd<-dplyr::select(predsm8c4Fd, lon, lat, ypred)%>%rename(x=lon, y=lat)
predsm9c4Fd<-dplyr::select(predsm9c4Fd, lon, lat, ypred)%>%rename(x=lon, y=lat)


#plot predictions----
#BEC plot locations
load(file="data/tree_data_cleaned.Rdata")
Hw<-subset(tree_dat, Species=='TSUGHET') 
Hw<-mutate(Hw, cover_rank=case_when(TotalA>20~1, TotalA>0 & TotalA<7.5~3,
                                    TRUE~2))
Hw1<-subset(Hw, cover_rank!="0")
Hw1zonal<-subset(Hw1, NutrientRegime_clean=="C"& MoistureRegime_clean=='4')
Hwzonal<-subset(Hw, NutrientRegime_clean=="C"& MoistureRegime_clean=='4')

Cw<-subset(tree_dat, Species=='THUJPLI') 
Cw<-mutate(Cw, cover_rank=case_when(TotalA>20~1, TotalA>0 & TotalA<7.5~3,
                                    TRUE~2))
Cw1<-subset(Cw, cover_rank=="1")

Fd<-subset(tree_dat, Species=='PSEUMEN') 
Fd<-mutate(Fd, cover_rank=case_when(TotalA>20~1, TotalA>0 & TotalA<7.5~3,
                                    TRUE~2))
Fd1<-subset(Fd, cover_rank=="1")

#load BC boundary for masking 
bcboundary<-bcmaps::bc_bound_hres()
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84 

#turn preds into raster & plot-m5---- 
predsm5c4Hwr<-rasterFromXYZ(predsm5c4Hw)
predsm5c4Hwr<-mask(predsm5c4Hwr, bcbound.reproj)#mask areas not in BC boundary 
predsm5c4Hwr$ypred <- as.factor(predsm5c4Hwr$ypred) #not working??
plot(predsm5c4Hwr)#still coded backward - but can interpret that green = primary 
plot(predsm5c4Hwr,main="Mod5 Hw c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Hw1$Longitude, Hw1$Latitude, #add BEC plots where Hw is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))

points(x = tree_dat_map$Longitude, tree_dat_map$Latitude, 
       col="lightblue",
       cex = 0.1)

#turn preds into raster & plot-m6----
predsm6c4Hwr<-rasterFromXYZ(predsm6c4Hw)
predsm6c4Hwr<-mask(predsm6c4Hwr, bcbound.reproj)#mask areas not in BC boundary 
predsm6c4Hwr$ypred <- as.factor(predsm6c4Hwr$ypred) #not working??
plot(predsm6c4Hwr)#still coded backward - but can interpret that green = primary 
plot(predsm6c4Hwr,main="Mod6 Hw c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Hw1$Longitude, Hw1$Latitude, #add BEC plots where Hw is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))

#turn preds into raster & plot-m7----
predsm7c4Hwr<-rasterFromXYZ(predsm7c4Hw)
predsm7c4Hwr<-mask(predsm7c4Hwr, bcbound.reproj)#mask areas not in BC boundary 
predsm7c4Hwr$ypred <- as.factor(predsm7c4Hwr$ypred) #not working??
plot(predsm7c4Hwr)#still coded backward - but can interpret that green = primary 

plot(predsm7c4Hwr,main="Mod7 Hw c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Hw1$Longitude, Hw1$Latitude, #add BEC plots where Hw is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))


#turn preds into raster & plot-m8----
#Hw
predsm8c4Hwr<-rasterFromXYZ(predsm8c4Hw)
predsm8c4Hwr<-mask(predsm8c4Hwr, bcbound.reproj)#mask areas not in BC boundary 
#predsm8c4Hwr$ypred <- as.factor(predsm8c4Hwr$ypred) #not working??
#plot(predsm8c4Hwr)#still coded backward - but can interpret that green = primary 

plot(predsm8c4Hwr,main="Hw c4 2km pred feasibility- 11 derived params", legend=)#turn off auto legend  
points(x = Hw1zonal$Longitude, Hw1zonal$Latitude, #add BEC plots where Hw is primary (area> 20%)
       col="black",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3" , "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "black"), 
       bty="n",cex=0.8, 
       inset=c(0.02, 0.05))

#Cw
predsm8c4Cwr<-rasterFromXYZ(predsm8c4Cw)
predsm8c4Cwr<-mask(predsm8c4Cwr, bcbound.reproj)#mask areas not in BC boundary 
#predsm8c4Cwr$ypred <- as.factor(predsm8c4Cwr$ypred) #not working??
#plot(predsm8c4Cwr)#still coded backward - but can interpret that green = primary 

plot(predsm8c4Cwr,main="Mod8 Cw c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Cw1$Longitude, Cw1$Latitude, #add BEC plots where Cw is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3"), #"BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1"),# "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))

#Fd
predsm8c4Fdr<-rasterFromXYZ(predsm8c4Fd)
predsm8c4Fdr<-mask(predsm8c4Fdr, bcbound.reproj)#mask areas not in BC boundary 
#predsm8c4Fdr$ypred <- as.factor(predsm8c4Fdr$ypred) #not working??
#plot(predsm8c4Fdr)#still coded backward - but can interpret that green = primary 

plot(predsm8c4Fdr,main="Mod8 Fd c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Fd1$Longitude, Fd1$Latitude, #add BEC plots where Fd is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))

#turn preds into raster & plot-m9----
#HW
predsm9c4Hwr<-rasterFromXYZ(predsm9c4Hw)
predsm9c4Hwr<-mask(predsm9c4Hwr, bcbound.reproj)#mask areas not in BC boundary 
#predsm9c4Hwr$ypred <- as.factor(predsm9c4Hwr$ypred) #not working??
#plot(predsm9c4Hwr)#still coded backward - but can interpret that green = primary 

plot(predsm9c4Hwr,main="Hw c4 2km pred feasibility-17 params", legend=F)#turn off auto legend  
points(x = Hw1zonal$Longitude, Hw1zonal$Latitude, #add BEC plots where Hw is primary (area> 20%)
       col="black",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "black"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))

#Cw
predsm9c4Cwr<-rasterFromXYZ(predsm9c4Cw)
predsm9c4Cwr<-mask(predsm9c4Cwr, bcbound.reproj)#mask areas not in BC boundary 
#predsm9c4Cwr$ypred <- as.factor(predsm9c4Cwr$ypred) #not working??
#plot(predsm9c4Cwr)#still coded backward - but can interpret that green = primary 

plot(predsm9c4Cwr,main="Mod9 Cw c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Cw1$Longitude, Cw1$Latitude, #add BEC plots where Cw is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))
#Fd
predsm9c4Fdr<-rasterFromXYZ(predsm9c4Fd)
predsm9c4Fdr<-mask(predsm9c4Fdr, bcbound.reproj)#mask areas not in BC boundary 
#predsm9c4Fdr$ypred <- as.factor(predsm9c4Fdr$ypred) #not working??
#plot(predsm9c4Fdr)#still coded backward - but can interpret that green = primary 

plot(predsm9c4Fdr,main="Mod9 Fd c4 2km predicted feasibility", legend=F)#turn off auto legend  
points(x = Fd1$Longitude, Fd1$Latitude, #add BEC plots where Fd is primary (area> 20%)
       col="lightgrey",
       cex = 0.1)
legend("topright",   
       legend =  c("1","2", "3", "BEC plot"), 
       # the legend
       fill =  c("green4","greenyellow", "orange1", "lightgrey"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))




#try with tidy terra
library(tidyterra)
predsm8c4Hwrx<-rast(predsm8c4Hwr)

ggplot() + 
    geom_spatraster(data = predsm8c4Hwrx) +
  #scale_fill_gradient(direction=-1)
scale_fill_grass_b(palette = "celsius", breaks = c(0, 1,2,3))+ ggtitle("Hw zonal 11 derived params")+
  geom_point(data=Hwzonal, aes(x= Longitude, y=Latitude), shape = 1,  alpha=0.5, color='red', size = .7)


#  scale_fill_gradient(colors = c("green4","greenyellow", "orange1", "lightgrey"), breaks = seq(0, 3, 1))+
#  scale_fill_gradient(palette = 'viridis', breaks = seq(0, 3, 1), alpha = 0.5)+

cols <- c("1" = "green4", "2" = "greenyellow", "3" = "orange1", "0" = "lightgrey")
  
scale_fill_discrete( )
scale_fill_gradient(
  colors = c("green4","greenyellow", "orange1", "lightgrey"))

  
  
  
  ,
  scale_colour_grass_d(
    palette = "viridis",
   discrete_scale(aesthetics = accent),
    alpha = 1,
    direction = 1,
    na.translate = FALSE,
    drop = TRUE
  )
  scale_fill_grass_d(
    direction = -1)

geom_spatraster(data = temp_rast, aes(fill = tavg_04))



