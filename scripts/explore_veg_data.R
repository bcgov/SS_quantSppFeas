#License info ----
#Copyright 2019 Province of British Columbia
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License  is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#libraries 
library(tidyverse)
require(data.table)

#pull in climate & plot data ----
load(file= 'data/clim_dat.plots.Rdata') #59314 obs

#merge with veg & site data
BEC_data<-readRDS("data/BEC_data.rds") #change this to OS #13 
veg_dat<-BEC_data$veg
all_dat<-left_join(veg_dat, plot_dat)#join with plot & climate info 

site_dat<-BEC_data$info
all_dat<-left_join(all_dat, rename(site_dat, PlotNumber=Plot))#join with plot & climate info 

#clean up veg data----
rm(plot_dat)
rm(site_dat)
rm(veg_dat)
gc()

#look at data across species/sites 
#match with updated site series info from Will - awaiting final BEC v13 x plot numbers list
ss_cleaned<-read.csv("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_workingfiles/Feasibility_modelling/All_BGC12DEC2024_SU.csv") 
all_dat<-left_join(all_dat, ss_cleaned) 

#remove anything that is not tree layer (TotalA or B)
#how to deal with B layer?? add to A layer? if nothing in A layer, use B layer? average? 
tree_dat<-filter(all_dat, !is.na(TotalA)|!is.na(TotalB))  
#(not designated A BGC Unit from Will's list) and #!is.na(bgc))
rm(all_dat)
gc()
tree_dat<-separate(tree_dat, SiteUnit, into = c("Zone", "Site"), sep = '/', remove=F) #can ignore warning 
sort(unique(tree_dat$Species))

spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>100)
spp_keep<-spp_keep$Species #30 total 


#remove any tree spp with <100 total obs  
tree_dat<-filter(tree_dat, Species %in% spp_keep)

#keep only major tree species 
tree_dat<-subset(tree_dat, grepl('ABIE|THUJ|PSEU|PINU|PICE|CALL|TAX|LARI|TSUG|ACER|POPU|ALNURUB|BETUPAP|
                                 ACERMAC|ARBU|QUER', Species))

#look at nobs of spp by site units 
spp_tab<-tree_dat%>%
                group_by(Species, SiteUnit, bgc,Site)%>%
                summarise(nobs=n())

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
tree_dat<-dplyr::select(tree_dat, all_of(cols_keep))

#remove info not using
tree_dat<-dplyr::select(tree_dat, 
    -FSRegionDistrict, -SV_FloodPlain, -ProvinceStateTerritory, 
    -NtsMapSheet, -Flag, -SpeciesListComplete, -UpdatedFromCards, -Zone, -Ecosection)
names(tree_dat)

#aspect, slope 
sort(names(tree_dat))
str(tree_dat)#what are continuous/numeric?? elevation, slope, aspect, substrate categories 
hist(tree_dat$Aspect)#Should be from 0-360
tree_dat<-subset(tree_dat, Aspect<361)
hist(tree_dat$SlopeGradient)#Should be from 0-100??
tree_dat<-subset(tree_dat, SlopeGradient<101)

# SNR, aSMR 
unique(tree_dat$NutrientRegime)# need to finalize calls on all
unique(tree_dat$MoistureRegime) # need to finalize calls on all 

#try to create table to fill in from
SMRs<-group_by(tree_dat, SiteUnit)%>%dplyr::select(MoistureRegime, NutrientRegime) %>%distinct()%>%mutate(ctMoist=n_distinct(MoistureRegime))%>%
  mutate(ctNut=n_distinct(NutrientRegime))%>%mutate(MoistNA= is.na(MoistureRegime))%>%mutate(NutNA= is.na(NutrientRegime))

#set rules for transitional SMRs and SNRs 
sort(table(SMRs$NutrientRegime)) #table-base R
                                
SMRs<-mutate(SMRs, NutrientRegime_clean = case_when(NutrientRegime=="A" ~ "A",
                                      NutrientRegime=="A+" ~ "A",
                                      NutrientRegime=="AB" ~ "B",
                                      NutrientRegime=="B" ~ "B",
                                      NutrientRegime=="BA" ~ "B",
                                      NutrientRegime=="B+" ~ "B",
                                      NutrientRegime=="B-" ~ "B",
                                      NutrientRegime=="BC" ~ "C",
                                      NutrientRegime=="C" ~ "C",
                                      NutrientRegime=="C+" ~ "C",
                                      NutrientRegime=="C-" ~ "C",
                                      NutrientRegime=="CB" ~ "C",
                                      NutrientRegime=="CD" ~ "C",
                                      NutrientRegime=="D" ~ "D",
                                      NutrientRegime=="D+" ~ "D",
                                      NutrientRegime=="D-" ~ "D",
                                      NutrientRegime=="DC" ~ "C",
                                      NutrientRegime=="DE" ~ "D",
                                      NutrientRegime=="E" ~ "E",
                                      NutrientRegime=="E-" ~ "E",
                                      NutrientRegime=="ED-" ~ "D",
                                      NutrientRegime=="F" ~ "F",
                                      NutrientRegime=="F-" ~ "F",
                                      NutrientRegime=="M" ~ "C",
                                      NutrientRegime=="P" ~ "B",
                                      NutrientRegime=="R" ~ "D",
                                      NutrientRegime=="3" ~ "NA",
                                      NutrientRegime=="6" ~ "NA"))

sort(table(SMRs$MoistureRegime))
SMRs<-mutate(SMRs, MoistureRegime_clean = case_when(MoistureRegime=="$" ~ "NA",
                                                    MoistureRegime=="0" ~ "0",
                                                    MoistureRegime=="0-1" ~ "1",
                                                    MoistureRegime=="0+" ~ "0",
                                                    MoistureRegime=="1" ~ "1",
                                                    MoistureRegime=="1-0" ~ "1",
                                                    MoistureRegime=="1-1" ~ "1",
                                                    MoistureRegime=="1-2" ~ "2",
                                                    MoistureRegime=="1+" ~ "1",
                                                    MoistureRegime=="2" ~ "2",
                                                    MoistureRegime=="2-" ~ "2",
                                                    MoistureRegime=="2-0" ~ "1",
                                                    MoistureRegime=="2-1" ~ "2",
                                                    MoistureRegime=="2-3" ~ "3",
                                                    MoistureRegime=="2+" ~ "2",
                                                    MoistureRegime=="3" ~ "3",
                                                    MoistureRegime=="3-" ~ "3",
                                                    MoistureRegime=="3-2" ~ "3",
                                                    MoistureRegime=="3-3" ~ "3",
                                                    MoistureRegime=="3-4" ~ "4",
                                                    MoistureRegime=="3(2" ~ "3",
                                                    MoistureRegime=="3+" ~ "3",
                                                    MoistureRegime=="4" ~ "4",
                                                    MoistureRegime=="4-" ~ "4",
                                                    MoistureRegime=="4-3" ~ "4",
                                                    MoistureRegime=="4-5" ~ "4",
                                                    MoistureRegime=="4(3" ~ "4",
                                                    MoistureRegime=="4.5" ~ "4", 
                                                    MoistureRegime=="4+" ~ "4",
                                                    MoistureRegime=="5" ~ "5",
                                                    MoistureRegime=="5-4" ~ "4",
                                                    MoistureRegime=="5-6" ~ "5",
                                                    MoistureRegime=="5(6" ~ "5",
                                                    MoistureRegime=="5H" ~ "5",
                                                    MoistureRegime=="6" ~ "6",
                                                    MoistureRegime=="6-" ~ "6",
                                                    MoistureRegime=="6-4" ~ "5",
                                                    MoistureRegime=="6-5" ~ "5",
                                                    MoistureRegime=="6-7" ~ "6", 
                                                    MoistureRegime=="6+" ~ "6",
                                                    MoistureRegime=="7" ~ "7",
                                                    MoistureRegime=="7-" ~ "7",
                                                    MoistureRegime=="8" ~ "8",
                                                    MoistureRegime=="PM" ~ "NA"))

SMRs<-dplyr::select(SMRs, SiteUnit, MoistureRegime, NutrientRegime, NutrientRegime_clean, MoistureRegime_clean)

#merge back to tree data
tree_dat<-left_join(tree_dat, SMRs)
unique(tree_dat$MoistureRegime_clean)
unique(tree_dat$NutrientRegime_clean)

#remove rows with NAs in Moisture/Nutrient regimes 
tree_dat<- subset(tree_dat, NutrientRegime_clean!="NA"& MoistureRegime_clean!="NA")

#look at other veg vars of interest
hist(tree_dat$StrataCoverTree) #normally distributed-ish  #is this all tree spp combined??
hist(tree_dat$StrataCoverTotal) #what is this?
hist(tree_dat$TotalA) #right skewed - beta dist?
hist(log(tree_dat$TotalA)) #not bad... 

#create a year column
#library(lubridate)
#tree_dat<-mutate(tree_dat,year=year(Date))
names(tree_dat)

#fix species names 
tree_dat<-mutate(tree_dat, Species=if_else(Species=="PINUCON1"|Species=="PINUCON2","PINUCON", Species))%>%
  mutate(Species=if_else(Species=="PSEUMEN1"|Species=="PSEUMEN2","PSEUMEN", Species))

unique(tree_dat$Species)

#check for NAs 
tree_dat2<-na.omit(tree_dat)
nas<-anti_join(tree_dat, tree_dat2)

nacols <- function(df) {
  colnames(df)[unlist(lapply(df, function(x) any(is.na(x))))]
}
nacols(nas)
names(tree_dat)

#save cleaned tree data----
save(tree_dat, file="data/tree_data_cleaned.Rdata")
rm(all_dat)
gc()

#add zeroes to plots where species not observed---- 
#load(file="data/tree_data_cleaned.Rdata")
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>250)
spp_keep<-spp_keep$Species #16 

#filter to 16
tree_dat<-subset(tree_dat, Species %in% spp_keep)
unique(tree_dat$Species)              

#expand grid
tree_dat_wzeros<-expand.grid(PlotNumber=unique(tree_dat$PlotNumber), Species=spp_keep)

#bring back in climate data by plot 
load(file= 'data/clim_dat.plots.Rdata')
tree_dat_wzeros<-left_join(tree_dat_wzeros, plot_dat)
gc()

#merge back in plot data (minus climate)
sort(names(tree_dat))
tree_dat<-dplyr::select(tree_dat,  PlotNumber, Species, TotalA, ID, ProjectID, Date, #SiteSurveyor,PlotRepresenting, Location,
  Longitude, Latitude, LocationAccuracy, SubZone, SiteSeries, MoistureRegime,           
  NutrientRegime, Elevation, SlopeGradient, Aspect, MesoSlopePosition, #SubstrateDecWood, SubstrateBedRock,
  #SubstrateRocks, SubstrateMineralSoil, SubstrateOrganicMatter, SubstrateWater, SurficialMaterialSurf, 
  SoilDrainage, HumusForm, StrataCoverTree, StrataCoverShrub, StrataCoverHerb, StrataCoverMoss, UserSiteUnit, GIS_BGC, GIS_BGC_VER,
  StrataCoverTotal, Elevation_overlay, SiteUnit, Site, bgc, NutrientRegime_clean, MoistureRegime_clean, year)
#pull out tree cover 
plot_dat2<-dplyr::select(tree_dat, -Species, -TotalA, -ID)%>%distinct(.)
tree_dat_wzeros<-left_join(tree_dat_wzeros, plot_dat2) #check warnings about duplicated info 
tree_dat_wzeros<-distinct(tree_dat_wzeros)
#put tree cover back in 
tree_dat2<-dplyr::select(tree_dat, PlotNumber, Species, TotalA, ID)%>%distinct(.)
tree_dat_wzeros<-left_join(tree_dat_wzeros, tree_dat2)%>%relocate(c(TotalA, ID), .after = Species)

#now add in zeroes where total A is NA
tree_dat_wzeros<-mutate(tree_dat_wzeros, TotalA=if_else(is.na(TotalA), 0, TotalA))

hist(tree_dat_wzeros$TotalA) #very zero inflated 

#check for NAs 
sort(nacols(tree_dat_wzeros))

#save
save(tree_dat_wzeros, file="data/tree_data_cleaned_wzeros.Rdata") #too big to push- save on OS #13


#look at most measured spp----- 
#Western red Cedar 
Cw<-subset(tree_dat, Species=='THUJPLI')  #3639 obs
ggplot(Cw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 
#geom_smooth(method='lm', se = F)#+
#theme(legend.position = 'none')
#xlim(-2.5, 12) #git rid of outliers 

#Western Hemlock
Hw<-subset(tree_dat, Species=='TSUGHET') #4716 obs
ggplot(Hw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)
#geom_smooth(method='lm', se = F)#+
#theme(legend.position = 'none')
# xlim(-2.5, 12) #git rid of outliers 

#Engelman Spruce
Se<-subset(tree_dat, Species=='PICEENE') #4055 obs
ggplot(Se, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Douglas Fir 
Fd<-subset(tree_dat, Species=="PSEUMEN") #3705 obs
ggplot(Fd, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Subalpine Fir
Bl<-subset(tree_dat, Species=="ABIELAS") #4060 obs 
ggplot(Bl, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Lodgepole Pine 
Pl<-subset(tree_dat, Species=="PINUCON") #3164 obs
ggplot(Pl, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Amabalis fir
Ba<-subset(tree_dat, Species=="ABIEAMA") #1913 obs
ggplot(Ba, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Mountain Hemlock
Hm<-subset(tree_dat, Species=="TSUGMER") #1337 obs 
ggplot(Hm, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+ 
  facet_wrap(~bgc) #why some blank here??

#Yellow Cedar
Yc<-subset(tree_dat, Species=="CALLNOO") #1187 obs
ggplot(Yc, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Sitka Spruce 
Ss<-subset(tree_dat, Species=="PICESIT") #939 obs (>1000 cut off?)
ggplot(Ss, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Trembling Aspen 
At<-subset(tree_dat, Species=="POPUTRE") #536 obs 
ggplot(At, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Paper birch 
Ep<-subset(tree_dat, Species=="BETUPAP") #531 obs 
ggplot(Ep, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#Ponderosa Pine
Py<-subset(tree_dat, Species=="PINUPON") #568 obs
ggplot(Py, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#White Spruce 
Sw<-subset(tree_dat, Species=="PICEGLA")  #397 obs (<500 cut off?)
ggplot(Sw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#Western Larch 
Lw<-subset(tree_dat, Species=="LARIOCC") #442 obs(<500 cut off?)
ggplot(Lw, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc) 

#Black Spruce 
Sb<-subset(tree_dat, Species=="PICEMAR") #313 obs (<500 cut off?)
ggplot(Sb, aes(y=TotalA, x=MAT, color=Site))+
  geom_point()+
  facet_wrap(~bgc)

#uneven sampling across BGCs by species, may need stronger threshold for number obs per BGC per site, species etc
#currently n=2+ 

#Save top 16 species level datasets
#save(Cw, Hw, Se, Fd, Bl, Pl, Ba, Hm, Yc, Ss, At, Ep, Py, Sw, Lw, Sb, file="data/tree_spp_data_cleaned.Rdata")



