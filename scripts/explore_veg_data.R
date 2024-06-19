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
source("scripts/climr_getdata.R") #ignore warnings

#merge with veg data---- 
veg_dat<-BEC_data$veg
all_dat<-left_join(veg_dat, plot_dat)#join with plot & climate info 

site_dat<-BEC_data$info
all_dat<-left_join(all_dat, rename(site_dat, PlotNumber=Plot))#join with plot & climate info 

#clean up veg data----
rm(plot_dat)
rm(site_dat)
rm(my_points)
rm(veg_dat)
gc()

#look at data across species/sites 
#match with updated site series info from Will 5/16/24
ss_cleaned<-read.csv('data/updated_siteseries.csv') #what defines why these plots are 'good'??
all_dat<-left_join(all_dat, ss_cleaned) 

#remove anything not designated A BEC or BGC Unit from Will's list and that is not tree layer (TotalA)
tree_dat<-filter(all_dat, !is.na(bgc) &!is.na(TotalA))
tree_dat<-separate(tree_dat, SiteUnit, into = c("Zone", "Site"), sep = '/', remove=F) #can ignore warning 
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

#look at nobs of spp by site units 
spp_tab<-tree_dat%>%
                group_by(Species, SiteUnit, bgc,Site)%>%
                summarise(nobs=n())
max(spp_tab$nobs)#228
min(spp_tab$nobs) #1
mean(spp_tab$nobs)#7
#filter anything with n obs <3 (need for random effect)
remove<-subset(spp_tab, nobs<2)
tree_dat<-anti_join(tree_dat, remove)
spp_tab<-anti_join(spp_tab, remove)

#look at enviro vars of interest
#which params are consistently measured?
col_cts<-as.data.frame(colSums(!is.na(tree_dat)))
#clean up
col_cts$cols<-row.names(col_cts)
col_cts$n_obs<-col_cts$`colSums(!is.na(tree_dat))`
col_cts$`colSums(!is.na(tree_dat))`<-NULL
cols_keep<-subset(col_cts, n_obs>30000)
cols_keep<-cols_keep$cols

#filter for only consistently recorded columns 
#tree_datx<-select(tree_dat, all_of(cols_keep))


#remove info not using
#tree_datx<-select(tree_datx, 
#    -Cover2, -FSRegionDistrict, -SV_FloodPlain, -ProvinceStateTerritory, 
#    -NtsMapSheet, -Flag, -SpeciesListComplete, -UpdatedFromCards, -Zone)
#names(tree_datx)

# SNR, aSMR, aspect, slope 
sort(names(tree_dat))
str(tree_dat)#what are continuous/numeric?? elevation, slope, aspect, substrate categories 
hist(tree_dat$Aspect)#Should be from 0-360
tree_dat<-subset(tree_dat, Aspect<361)
hist(tree_dat$SlopeGradient)#Should be from 0-100??
tree_dat<-subset(tree_dat, SlopeGradient<101)


#look at other veg vars of interest
hist(tree_dat$StrataCoverHerb)
hist(tree_dat$StrataCoverMoss)
hist(tree_dat$StrataCoverShrub)
hist(tree_dat$StrataCoverTree) #normally distributed-ish  #is this all tree spp combined??
hist(tree_dat$StrataCoverTotal) #what is this?


hist(tree_dat$TotalA) #right skewed - beta dist?

#create a year column
library(lubridate)
tree_dat<-mutate(tree_dat,year=year(Date))


#save(tree_dat, file="data/tree_data_cleaned.Rdata")

#read in cleaned tree data---- 
load(file="data/tree_data_cleaned.Rdata")

#bring in feas tables---- 
feas_tab<-read.csv("data/FeasibilityUpdates.csv")#downloaded from ByBEC 6/3/24
feas_tab2<-read.csv("data/Feasibility_v12_15.csv")#most recent version from ccisr/data-raw/data_tables/

#check that feasibilities same between BYBEC and ccissr versions----
names(feas_tab)
feas_tab2<-select(feas_tab2, bgc, ss_nospace, sppsplit, feasible, newfeas, mod)%>%rename(spp=sppsplit)
diff<-anti_join(feas_tab, feas_tab2)
diff2<-anti_join(feas_tab2, feas_tab)

feas_tab<-anti_join(feas_tab, diff)
feas_tab<-rbind(feas_tab, diff2)

feas_tab2<-rename(feas_tab2, newfeas2=newfeas)
allfeas<-full_join(feas_tab, feas_tab2)
cor(allfeas$newfeas, allfeas$newfeas2) 

allfeas<-mutate(allfeas, check= if_else(newfeas2==newfeas, T, F)) #only 2 different
#CWHms1 CWHms1/03 Ba 3 3 SAS-HAK 2 FALSE
#CWHms1 CWHms1/03 Ba 3 2 SAS-HAK 3 FALSE



#update species naming----
unique(sort(feas_tab2$spp))
unique(sort(feas_tab$spp))

feas_tab2<-mutate(feas_tab2, Species= case_when(spp=="Ba"~"ABIEAMA", 
                                              spp=="Bg"~"ABIEGRA", 
                                              spp=="Bl"~"ABIELAS", 
                                              spp=="Mb"~"ACERMAC", 
                                              spp=="Dr"~"ALNURUB",
                                              spp=="Ra"~"ARBUMEN",  
                                              spp=="Ep"~"BETUPAP", 
                                              spp=="Yc"~"CALLNOO", 
                                              spp=="Lw"~"LARIOCC",  
                                              spp=="Se"~"PICEENE",
                                              spp=="Sw"~"PICEGLA",
                                              spp=="Sb"~"PICEMAR",
                                              spp=="Ss"~"PICESIT",
                                              spp=="Sxl"~"PICEXLU",
                                              spp=="Pa"~"PINUALB",
                                              spp=="Pl"~"PINUCON",
                                              spp=="Plc"~"PINUCON", #coastal
                                              spp=="Pli"~"PINUCON", #interior
                                              spp=="Pyc"~"PINUPON",
                                              spp=="Pyi"~"PINUPON",
                                              spp=="Acb"~"POPUBAL",
                                              spp=="At"~"POPUTRE",
                                              spp=="Act"~"POPUTRI",
                                              spp=="Fd"~"PSEUMEN", 
                                              spp=="Fdi"~"PSEUMEN", #interior
                                              spp=="Fdc"~"PSEUMEN", #coastal
                                              spp=="Tw"~ "TAXUBRE",
                                              spp=="Cw"~ "THUJPLI",
                                              spp=="Hw"~ "TSUGHET",
                                              spp=="Hm"~"TSUGMER")) 


tree_dat<-mutate(tree_dat, Species=if_else(Species=="PINUCON1"|Species=="PINUCON2","PINUCON", Species))%>%
  mutate(Species=if_else(Species=="PSEUMEN1"|Species=="PSEUMEN2","PSEUMEN", Species))

sort(unique(feas_tab2$Species))
sort(unique(tree_dat$Species))
#have plot data but missing feasibility ratings on Pw (PINUMON-Western white pine), Tw (TAXUBRE-Western Yew) 

#combine feas table with plot data----
tree_datx<-mutate(tree_dat,  ss_nospace= gsub(" ", "", SiteUnit)) #create matching column to feas table
tree_datx<-left_join(tree_datx, feas_tab2, relationship = "many-to-many")%>%rename(newfeas=newfeas2)
#need to fix the 01 to 101s in several places so these can align correctly 

#look at whether feasibility is reflective of plot level abundance by species 
ggplot(tree_datx, aes(y=TotalA, x=newfeas))+
  geom_point()+
  facet_wrap(~Species) 

cor.test(tree_datx$TotalA, tree_datx$newfeas) #-0.25

cors<- group_by(tree_datx, Species)%>% 
  summarise(abun_feas_cor=cor(TotalA, newfeas, use="na.or.complete")) #-0.75 to +0.25 varies a lot by spp 
#two spp with pos correlations which is opposite of expected-> poor data coverage  
#POPUTRI 0.24638427
#PINUALB 0.14793807 

#look at most measured spp- 
#Western red Cedar 
Cw<-subset(tree_datx, Species=='THUJPLI')  #3814 obs
ggplot(Cw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 
#geom_smooth(method='lm', se = F)#+
  #theme(legend.position = 'none')
  #xlim(-2.5, 12) #git rid of outliers 

#Western Hemlock
Hw<-subset(tree_datx, Species=='TSUGHET') #4716 obs
ggplot(Hw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)
  #geom_smooth(method='lm', se = F)#+
  #theme(legend.position = 'none')
 # xlim(-2.5, 12) #git rid of outliers 

#Engelman Spruce
Se<-subset(tree_datx, Species=='PICEENE') #4055 obs
ggplot(Se, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Douglas Fir 
Fd<-subset(tree_datx, Species=="PSEUMEN") #3705 obs
ggplot(Fd, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Subalpine Fir
Bl<-subset(tree_datx, Species=="ABIELAS") #4060 obs 
ggplot(Bl, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Lodgepole Pine 
Pl<-subset(tree_datx, Species=="PINUCON") #3164 obs
ggplot(Pl, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Amabalis fir
Ba<-subset(tree_datx, Species=="ABIEAMA") #1913 obs
ggplot(Ba, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Mountain Hemlock
Hm<-subset(tree_datx, Species=="TSUGMER") #1337 obs 
ggplot(Hm, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+ 
  facet_wrap(~bgc) #why some blank here??

#Yellow Cedar
Yc<-subset(tree_datx, Species=="CALLNOO") #1187 obs
ggplot(Yc, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Sitka Spruce 
Ss<-subset(tree_datx, Species=="PICESIT") #939 obs (>1000 cut off?)
ggplot(Ss, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Trembling Aspen 
At<-subset(tree_datx, Species=="POPUTRE") #536 obs 
ggplot(At, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Paper birch 
Ep<-subset(tree_datx, Species=="BETUPAP") #531 obs 
ggplot(Ep, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Ponderosa Pine
Py<-subset(tree_datx, Species=="PINUPON") #568 obs
ggplot(Py, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#White Spruce 
Sw<-subset(tree_datx, Species=="PICEGLA")  #397 obs (<500 cut off?)
ggplot(Sw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#Western Larch 
Lw<-subset(tree_datx, Species=="LARIOCC") #442 obs(<500 cut off?)
ggplot(Lw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#Black Spruce 
Sb<-subset(tree_datx, Species=="PICEMAR") #313 obs (<500 cut off?)
ggplot(Sb, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#uneven sampling across BGCs by species, may need stronger threshold for number obs per BGC per site, species etc
#currently n=2+ 
