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

#look at data across species/sites 
#match with updated site series info from Will 5/16/24
ss_cleaned<-read.csv('data/updated_siteseries.csv')
all_dat<-left_join(all_dat, ss_cleaned) 

#remove anything not designated A BEC or BGC Unit from Will's list and that is not tree layer (TotalA)
tree_dat<-filter(all_dat, !is.na(bgc) &!is.na(TotalA))
tree_dat<-separate(tree_dat, SiteUnit, into = c("Zone", "Site"), sep = '/', remove=F) #ignore warning 
sort(unique(tree_dat$Species))

#remove anything not assigned to spp level 
#tree_dat<-subset(tree_dat, Species!="SALIX" & Species!="TAXUS"&Species!="MALUS"&Species!="UNKNOWN"&Species!="PICEA"&
#     Species!="POPULUS"& Species!="BETULA" & Species!="ALNUS"  & Species!="TSUGA" & Species!="THUJA")
  
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>100)
spp_keep<-spp_keep$Species #30 total 

#remove any tree spp with <100 total obs  
tree_dat<-filter(tree_dat, Species %in% spp_keep) 

spp_tab<-tree_dat%>%
                group_by(Species,SiteUnit, bgc,Site)%>%
                summarise(nobs=n())

#look at most measured spp- 
#Western red Cedar 
Cw<-subset(tree_dat, Species=='THUJPLI')
ggplot(Cw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) #+
  #geom_smooth(method='lm', se = F)#+
  #theme(legend.position = 'none')
  #xlim(-2.5, 12) #git rid of outliers 

#Western Hemlock
Hw<-subset(tree_dat, Species=='TSUGHET') 
ggplot(Hw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)
  #geom_smooth(method='lm', se = F)#+
  #theme(legend.position = 'none')
 # xlim(-2.5, 12) #git rid of outliers 


#uneven sampling across BGCs, need some threshold for number of siteseries per BGC etc, also probably number of records per Site Series
all_dat<-filter(all_dat, !is.na(BECSiteUnit))


#add feasibility tables in??

names(all_dat)

#test model
library(lme4)
library(lmerTest)
library(dplyr)
hist(Hw$TotalA)

hist(log(Hw$TotalA))
hist(Hw$MAT)

unique(Hw$MAT)

#categorize MAT            
Hw<-mutate(Hw, MATcat=cut(MAT, breaks = 5,
      labels = c("MAT1", "MAT2", "MAT3", "MAT4", "MAT5")))


Hwmod<-lmer(log(TotalA)~ scale(MAT) + scale(PPT) + (MATcat|bgc:Site),  Hw) #can't do nested Random slope 
summary(Hwmod)

Hwmod<-lmer(log(TotalA)~ scale(MAT) + scale(PPT) + (MATcat|GIS_BGC) + (MATcat|GIS_BGC:SiteSeries),  Hw) #can't do nested RE

#climate, SNR, aSMR, aspect, slope 
