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


#libraries----
library(climr)
library(tidyverse)
library(terra)
library(data.table)

## provide a data.frame or data.table of point coordinates, IDs and elevation
BEC_data<-readRDS("data/BEC_data.rds")

#pull out plot data 
plot_dat<-BEC_data$env #70,547 plots

#make dataframe for extracting climate data
my_points <- select(plot_dat, Longitude, Latitude, Elevation, PlotNumber, ProjectID) %>%
  rename(lon = Longitude,   lat = Latitude, 
  elev = Elevation, id = PlotNumber)%>%
  na.omit() #remove NAs

#look at options 
#what to select here?
list_historic()
list_historic_ts()
list_variables() 
data('variables')#look up table for vars 

## climr query for the historic data - only using 1961-1990 for now 
clim_dat <- climr_downscale(
  xyz = my_points, which_normal = "auto",
  #historic_period = "2001_2020", 
  #historic_ts = C(1961:1990),
  #gcm_models = c("GFDL-ESM4", "EC-Earth3"), # specify two global climate models
  #ssp = c("ssp370", "ssp245"), # specify two greenhouse gas concentration scenarios
  #gcm_period = c("2001_2020", "2041_2060"), # specify two 20-year periods
  max_run = 3, # specify 3 individual runs for each model
  vars = c("PPT", "MAT")
)


#sanity check
#merge back with plot data 
plot_dat<-left_join(plot_dat, rename(clim_dat, PlotNumber=id))

ggplot(plot_dat, aes(x=Elevation, y=MAT))+
  geom_point()+
  geom_smooth(method='lm')+
  theme(legend.position = 'none')+
  ylim(-10, 10) + xlim(0,4000) #git rid of outliers 
#looks good, MAT does in fact decline with elevation

ggplot(plot_dat, aes(x=Elevation, y=PPT))+
  geom_point()+
  geom_smooth()+
  #ylim(-10, 10) + 
  xlim(0,4000) #git rid of outliers 
#precip more interesting - is this picking up snow??


#rasterize----
#using example from vignette https://bcgov.github.io/climr/articles/climr_workflow_beg.html#working-with-raster-data
## get the sample digital elevation model (dem) provided with `climr`
#dem_vancouver <- get(data("dem_vancouver")) |> 
#  unwrap()

## convert the DEM to a data.frame
#my_grid <- as.data.frame(dem_vancouver, cells = TRUE, xy = TRUE)
#colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

## A simple climr query. This will return the observed 1961-1990 and 2001-2020 mean annual temperature (MAT) for the raster grid points. 
#ds_out <- climr_downscale(
#  xyz = my_grid, 
#  historic_period = "2001_2020", 
#  vars = c("MAT")
#)

## populate the raster grid with the downscaled climate values
#my_clim <- rast(dem_vancouver) # use the DEM as a template raster
#my_clim[ds_out[PERIOD == "2001_2020", id]] <- ds_out[PERIOD == "2001_2020", MAT] # populate the raster cells with the 2001-2020 mean annual temperature (MAT) values, using the `id` field as the link. 

#plot(my_clim, main = "2001-2020 mean annual temperature (MAT)")
