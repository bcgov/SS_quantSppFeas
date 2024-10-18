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
BEC_data<-readRDS("data/BEC_data.rds") #change this to OS #13 

#pull out plot data 
plot_dat<-BEC_data$env #70,547 plots

#make dataframe for extracting climate data
my_points <- select(plot_dat, Longitude, Latitude, Elevation, PlotNumber, ProjectID) %>%
  rename(lon = Longitude,   lat = Latitude, 
  elev = Elevation, id = PlotNumber)%>%
  na.omit() #remove NAs

#look at options 
#what to select here?
list_obs_periods()
list_obs_years()
list_vars() 

vars<-climr::variables #look up table for vars 
var_names<-vars$Code

## climr query for the historic data - only using 1961-1990 for now 
## what is the resolution/scale of these data? PRISM 800m downscaled to plot-level (accuracy of GPS points and elevation- double checks elev vals make bigger difference)
cache_clear()
clim_dat <- downscale(
  xyz = my_points, which_refmap = "refmap_climr", 
  #historic_period = "2001_2020", 
  #historic_ts = C(1961:1990),
  #gcm_models = c("GFDL-ESM4", "EC-Earth3"), # specify two global climate models
  #ssp = c("ssp370", "ssp245"), # specify two greenhouse gas concentration scenarios
  #gcm_period = c("2001_2020", "2041_2060"), # specify two 20-year periods
  #max_run = 3, # specify 3 individual runs for each model
#  vars = c("PPT", "MAT", "CMD", 'AHM', 'CMI', 'DD5', 'TD', "PPT_10"))  #TD variable?? continentality??
 vars=var_names) 

#decide any other climate variables we want to include! 
#add derived variables used in BGC projections 
source("scripts/addVars.R")
addVars(clim_dat)

#sanity check
#merge back with plot data 
plot_dat<-left_join(plot_dat, rename(clim_dat, PlotNumber=id))
fwrite(plot_dat, "data/plot_dat_climr.csv")
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

ggplot(plot_dat, aes(x=Elevation, y=CMD))+
  geom_point()+
  geom_smooth()+
  #ylim(-10, 10) + 
  xlim(0,4000) #git rid of outliers 
#precip more interesting - is this picking up snow??
