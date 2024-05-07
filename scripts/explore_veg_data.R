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

#pull in climate data ----
source("scripts/climr_getdata.R")

#merge with veg data---- 
veg_dat<-BEC_data$veg
veg_dat<-left_join(veg_dat, rename(clim_dat, PlotNumber=id))#ignore warning 
names(veg_dat)
str(veg_dat)#what are the different covers (1-7) ?? 
unique(veg_dat$Species)

#sanity check
#look at cedar - should have pos relationship with PPT??  
Cw<-subset(veg_dat, Species=='THUJPLI')
ggplot(Cw, aes(y=Cover5, x=PPT))+
  geom_point()+
  geom_smooth(method='lm')#+
  #ylim(-2.5, 11) #git rid of outliers 
