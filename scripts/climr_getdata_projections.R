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

#Use PRISM DEM to call in downscaled climr data----
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="") #must have VPN connected- 800 m DEM
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))

## convert the DEM to a data.frame
my_grid <- as.data.frame(dem.bc, cells = TRUE, xy = TRUE)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

## climr call- This will return the observed 1961-1990 climates for the raster grid points. 
#which variables do we want? 
varsl = c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", 
         "CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp", 
         "PPT_sp", "PPT_sm", "Tmax_at", "Tmax_wt", "Tmax_sp", "Tmax_sm","Tmin_at", "Tmin_wt", "Tmin_sp", "Tmin_sm")
varsl<-unique(varsl)#make sure no dups 

cache_clear()
clim.bcv <- downscale(
  xyz = my_grid,  out_spatial = TRUE,   which_refmap = "refmap_climr",
  #obs_periods = "2001_2020", 
  vars = varsl)

saveRDS(clim.bcv, file="data/spatial/clim.bcv.RData")
plot(clim.bcv)
gc()

#rasterize from vector 
#clim.bcv<-readRDS(file="data/spatial/clim.bcv.RData")
my_rast <- rast(dem.bc) # use the DEM as a template raster
vars = c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", 
         "CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp", "PPT_sp", "PPT_sm", 'Tmax_at', "Tmax_wt", "Tmax_sp", "Tmax_sm",
         'Tmin_at', "Tmin_wt", "Tmin_sp", "Tmin_sm") 
clim.bcr<-rasterize(clim.bcv, my_rast, field= vars) #set background=0 for all cells that are NAs 
saveRDS(clim.bcr, file="data/spatial/clim.bcr.RData")
#clim.bcr<-readRDS(file="data/spatial/clim.bcr.RData")
plot(clim.bcr$TD)#looks good 
gc()

#where are the NAs? in the Ocean
#clim.bc <- as.data.frame(clim.bcr, cells = TRUE, xy = TRUE)
#clim.bcwna <- as.data.frame(clim.bcr, cells = TRUE, xy = TRUE, na.rm=F)
#climnas<-anti_join(clim.bcwna, clim.bc) #129284 NAs 
#climnas<-rename(climnas, id=cell, lon=x, lat=y)%>%dplyr::select(id, lat, lon)
#plot(clim.bcr$DD5)
#points(x = climnas$lon, 
#        y = climnas$lat, 
#        col = "red", 
#        cex = 0.5)
#write.csv(climnas, 'data/spatial/gpsnas.csv')

#aggregate up to 2km resolution
clim.bcr2k <- aggregate(clim.bcr, fact=3, fun=mean)
plot(clim.bcr2k)#looks good - can't really see resolution change
gc()

#clip to BC boundary
bcboundary<-bcmaps::bc_bound_hres()
#st_crs(bcboundary)
plot(st_geometry(bcboundary))
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84  
plot(st_geometry(bcbound.reproj))
clim.bcr2kcrop<-terra::crop(x = clim.bcr2k, y = bcbound.reproj) #crop climr output to BC boundary 
#saveRDS(clim.bcr2kcrop, file="data/spatial/clim.bcr2k.RData")
#clim.bcr2kcrop<-readRDS(file="data/spatial/clim.bcr2k.RData")
plot(clim.bcr2kcrop$DD5)#looks good- cropped slightly
gc()


#pull back out to df
clim.bc <- as.data.frame(clim.bcr2kcrop, cells = TRUE, xy = TRUE, na.rm=F)
clim.bc<-rename(clim.bc, id=cell, lon=x, lat=y)

#additional vars needed -"DD_delayed", "PPT_MJ", "PPT_JAS", "CMD.total", "CMDMax"----
clim.bc<-as.data.table(clim.bc) #requires data table 
source("scripts/addVars.R") 
addVars(clim.bc) #why not showing additional variables in environment??
save(clim.bc, file="data/clim.bc2k.RData")
