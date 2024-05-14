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
source("scripts/climr_getdata.R") #ignore warning

#merge with veg data---- 
veg_dat<-BEC_data$veg
all_dat<-left_join(veg_dat, plot_dat)#join with plot & climate info 

site_dat<-BEC_data$info
all_dat<-left_join(all_dat, rename(site_dat, PlotNumber=Plot))#join with plot & climate info 

#clean up
rm(plot_dat)
rm(site_dat)
rm(my_points)
rm(veg_dat)
gc()

sort(names(all_dat))
unique(all_dat$BECSiteUnit)

#look at data across species/sites 
all_dat<-filter(all_dat, !is.na(BECSiteUnit))

spp_tab<-all_dat%>%
                group_by(Species, BECSiteUnit, GIS_BGC, SiteSeries, UserSiteUnit)%>%
                summarise(nobs=n())

#subset to most measured spp- Western Hemlock
Hw<-subset(all_dat, Species=='TSUGHET') 
#plot across different sites 
#Site series are the different sites wet to dry (101-Zonal)- Needs cleaning/renumbering 
str(Hw)

ggplot(Hw, aes(y=TotalA, x=MAT, color=SiteSeries))+
  geom_point()+
  facet_wrap(~GIS_BGC)+
  geom_smooth(method='lm', se = F)+
  #theme(legend.position = 'none')
  xlim(-2.5, 12) #git rid of outliers 

#Doug Fir
Fd<-subset(all_dat, Species=='PSEUMEN1')
ggplot(Fd, aes(y=TotalA, x=MAT, color=SiteSeries))+
  geom_point()+
  facet_wrap(~GIS_BGC)+
  geom_smooth(method='lm', se = F)+
  #theme(legend.position = 'none')
  xlim(-2.5, 12) #git rid of outliers 

#Western red Cedar 
Cw<-subset(all_dat, Species=='THUJPLI')
ggplot(Cw, aes(y=TotalA, x=MAT, color=SiteSeries))+
  geom_point()+
  facet_wrap(~GIS_BGC)+
  geom_smooth(method='lm', se = F)+
  #theme(legend.position = 'none')
  xlim(-2.5, 12) #git rid of outliers 

#very uneven sampling across BGCs, need some threshold for number of siteseries per BGC etc, also probably number of records per Site Series
