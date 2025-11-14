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
#require(data.table)

#pull in climate & plot data ----
load(file= 'data/clim_dat.plots.Rdata') #59314 obs
#plot_dat<-plot_dat[ ,c(1:38)]

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

#look at data across site series 
#match with updated site series info from Will - awaiting final BEC v13 x plot numbers list
ss_cleaned<-read.csv("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_workingfiles/Feasibility_modelling/All_BGC12DEC2024_SU.csv") 
ss_cleaned<-read.csv("data/All_BGC13_May2025_SU.csv")
all_dat<-left_join(all_dat, ss_cleaned) 

#remove anything that is not tree layer (TotalA or B)
#how to deal with B layer?? add to A layer? if nothing in A layer, use B layer? average? 
tree_dat<-filter(all_dat, !is.na(TotalA)|!is.na(TotalB))  
rm(all_dat)
gc()
sort(unique(tree_dat$Species))

#keep only major tree species 
tree_dat<-subset(tree_dat, 
      grepl('THUJPLI|TSUGHET|PICEEN|PSEUMEN|ABIELAS|PINUCON|ABIEAMA|TSUGMER|CALLNOO|PICESIT|POPUTRE|POPUTRI|POPUBAL|
            BETUPAP|PINUPON|PICEGLA|LARIOCC|PICEMAR|PINUMON|ABIEGRA', Species))
#fix species names 
tree_dat<-mutate(tree_dat, Species=if_else(Species=="PINUCON1"|Species=="PINUCON2","PINUCON", Species))%>%
  mutate(Species=if_else(Species=="PSEUMEN1"|Species=="PSEUMEN2","PSEUMEN", Species))
tree_dat_sub<-mutate(tree_dat_sub, spp= case_when(Species=="ABIEAMA" ~ "Ba",
                                                  Species== "ARBUMEN" ~"Ra",
                                                  Species=="ABIEGRA"~ "Bg", 
                                                  Species=="ABIELAS"~"Bl", 
                                                  Species=="ACERMAC"~"Mb",
                                                  Species=="ACERCIR"~"Mv",
                                                  Species=="ALNURUB"~"Dr",
                                                  Species=="BETUPAP"~"Ep", 
                                                  Species=="CALLNOO"~"Yc",
                                                  Species=="LARILYA"~"La",
                                                  Species=="LARILAR"~"Lt",  
                                                  Species=="LARIOCC"~"Lw",  
                                                  Species=="PICEENE"~ "Sx", #not differentiated
                                                  Species=="PICEGLA"~"Sx", #not differentiated
                                                  Species=="PICEA"~"Sx", #not differentiated
                                                  Species=="PICEAX"~"Sx", #not differentiated
                                                  Species=="PICEENG"~"Sx", #not differentiated
                                                  Species=="PICEMAR"~"Sb",
                                                  Species=="PICESIT"~"Ss",#not differentiated
                                                  Species=="PICEXLU"~"Ss",#not differentiated
                                                  Species=="PINUALB"~"Pa",
                                                  Species=="PINUCON"~"Pl",#not differentiated
                                                  Species=="PINUPON"~"Py",
                                                  Species=="PINUMON"~"Pw",
                                                  Species=="POPUTRE" ~"At",
                                                  Species=="POPUBAL"~"Ac", #not differentiated
                                                  Species=="POPUTRI"~"Ac", #not differentiated
                                                  Species=="PSEUMEN"~"Fd", #not differentiated 
                                                  Species=="QUERGAR"~"Qg",
                                                  Species=="TAXUBRE" ~"Tw",
                                                  Species=="THUJPLI"~"Cw",
                                                  Species=="TSUGHET" ~"Hw",
                                                  Species=="TSUGMER"~"Hm")) 

#how many obs per spp?
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_tab0<-subset(spp_tab0, nobs>1000) #remove unclear hybrids
spp_keep<-spp_tab0$Species

tree_dat<-filter(tree_dat, Species %in% spp_keep)

#look at nobs of spp by site units 
spp_tab<-tree_dat%>%
                group_by(Species, SiteUnit)%>%
                summarise(nobs=n())

#filter anything with n obs <3 (need for random effect)
#remove<-subset(spp_tab, nobs<2)
#tree_dat<-anti_join(tree_dat, remove)
#spp_tab<-anti_join(spp_tab, remove)

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
#tree_dat<-dplyr::select(tree_dat, 
#    -FSRegionDistrict, -SV_FloodPlain, -ProvinceStateTerritory, 
#    -NtsMapSheet, -Flag, -SpeciesListComplete, -UpdatedFromCards, -Zone, -Ecosection)
#names(tree_dat)

#aspect, slope 
#sort(names(tree_dat))
#str(tree_dat)#what are continuous/numeric?? elevation, slope, aspect, substrate categories 
#hist(tree_dat$Aspect)#Should be from 0-360
#tree_dat<-subset(tree_dat, Aspect<361)
#hist(tree_dat$SlopeGradient)#Should be from 0-100??
#tree_dat<-subset(tree_dat, SlopeGradient<101)

# SNR, aSMR 
unique(tree_dat$NutrientRegime)# need to finalize calls on all
unique(tree_dat$MoistureRegime) # need to finalize calls on all 


#set rules for transitional SMRs and SNRs 
tree_dat<-mutate(tree_dat, NutrientRegime_clean = case_when(NutrientRegime=="A" ~ "A",
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

tree_dat<-mutate(tree_dat, MoistureRegime_clean = case_when(MoistureRegime=="$" ~ "NA",
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

#merge back to tree data
unique(tree_dat$MoistureRegime_clean)
unique(tree_dat$NutrientRegime_clean)

#remove rows with NAs in Moisture/Nutrient regimes 
tree_dat<- subset(tree_dat, NutrientRegime_clean!="NA"& MoistureRegime_clean!="NA") #2k rows 

#filter out poor quality 
tree_dat<-subset(tree_dat, !grepl("poor|Poor|POOR", SitePlotQuality))
tree_dat<-subset(tree_dat, !grepl("omit|Omit", UserSiteUnit))

#look at other veg vars of interest
hist(tree_dat$StrataCoverTree) #normally distributed-ish  #is this all tree spp combined??
hist(tree_dat$StrataCoverTotal) #what is this?
hist(tree_dat$TotalA) #right skewed - beta dist?
hist(log(tree_dat$TotalA)) #not bad... 

#create a year column
library(lubridate)
tree_dat<-mutate(tree_dat,year=year(Date))
names(tree_dat)

unique(tree_dat$Species)

#save cleaned tree data----
save(tree_dat, file="data/tree_data_cleaned.Rdata")
rm(BEC_data)
rm(clim_dat)
gc()

#add zeroes to plots where species not observed---- 
#NOT CURRENTLY USING 6/2025
#load(file="data/tree_data_cleaned.Rdata")
#spp_tab0<-tree_dat%>%
#  group_by(Species)%>%
#  summarise(nobs=n())
#spp_keep<-subset(spp_tab0, nobs>250)
#@spp_keep<-spp_keep$Species  

#filter 
#tree_dat<-subset(tree_dat, Species %in% spp_keep)
#unique(tree_dat$Species)              

#expand grid
#tree_dat_wzeros<-expand.grid(PlotNumber=unique(tree_dat$PlotNumber), Species=spp_keep)

#bring back in climate data by plot 
#load(file= 'data/clim_dat.plots.Rdata')
#tree_dat_wzeros<-left_join(tree_dat_wzeros, plot_dat)
#gc()

#merge back in plot data (minus climate)
#sort(names(tree_dat))
#tree_dat<-dplyr::select(tree_dat,  PlotNumber, Species, TotalA, TotalB, ID, ProjectID, Date, #SiteSurveyor,PlotRepresenting, Location,
 # Longitude, Latitude, LocationAccuracy, SubZone, SiteSeries, MoistureRegime,           
  #NutrientRegime, Elevation, SlopeGradient, Aspect, MesoSlopePosition, #SubstrateDecWood, SubstrateBedRock,
  #SubstrateRocks, SubstrateMineralSoil, SubstrateOrganicMatter, SubstrateWater, SurficialMaterialSurf, 
  #SoilDrainage, HumusForm, StrataCoverTree, StrataCoverShrub, StrataCoverHerb, StrataCoverMoss, UserSiteUnit, GIS_BGC, GIS_BGC_VER,
  #StrataCoverTotal, Elevation_overlay, NutrientRegime_clean, MoistureRegime_clean, year)

#pull out tree cover 
#plot_dat2<-dplyr::select(tree_dat, -Species, -TotalA, -TotalB, -ID)%>%distinct(.)
#tree_dat_wzeros<-left_join(tree_dat_wzeros, plot_dat2) #check warnings about duplicated info 
#tree_dat_wzeros<-distinct(tree_dat_wzeros)

#put tree cover back in 
#tree_dat2<-dplyr::select(tree_dat, PlotNumber, Species, TotalA, TotalB, ID)%>%distinct(.)
#tree_dat_wzeros<-left_join(tree_dat_wzeros, tree_dat2)%>%relocate(c(TotalA, ID), .after = Species)

#gc()

#now add in zeroes where cover is NA
#tree_dat_wzeros<-mutate(tree_dat_wzeros, TotalA=if_else(is.na(TotalA), 0, TotalA))%>%mutate(TotalB=if_else(is.na(TotalB), 0, TotalB))

#hist(tree_dat_wzeros$TotalA) #very zero inflated 
#hist(tree_dat_wzeros$TotalB) #very zero inflated 

#save
#save(tree_dat_wzeros, file="data/tree_data_cleaned_wzeros.Rdata") #too big to push- save on OS #13


#plot most measured spp----- 
load(file="data/tree_data_cleaned.Rdata") #only use actual plot data, not imputed zeros 

#Western red Cedar 
Cw<-subset(tree_dat, Species=='THUJPLI')  
ggplot(Cw, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC) 

#Western Hemlock
Hw<-subset(tree_dat, Species=='TSUGHET') 
ggplot(Hw, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC) 

#Engelman Spruce
Se<-subset(tree_dat, Species=='PICEENE')
ggplot(Se, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)


#Douglas Fir 
Fd<-subset(tree_dat, Species=="PSEUMEN") 
ggplot(Fd, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Subalpine Fir
Bl<-subset(tree_dat, Species=="ABIELAS") 
ggplot(Bl, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Lodgepole Pine 
Pl<-subset(tree_dat, Species=="PINUCON") 
ggplot(Pl, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Amabalis fir
Ba<-subset(tree_dat, Species=="ABIEAMA") 
ggplot(Ba, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Mountain Hemlock
Hm<-subset(tree_dat, Species=="TSUGMER")
ggplot(Hm, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+ 
  facet_wrap(~GIS_BGC) #why some blank here??

#Yellow Cedar
Yc<-subset(tree_dat, Species=="CALLNOO") 
ggplot(Yc, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Sitka Spruce 
Ss<-subset(tree_dat, Species=="PICESIT")
ggplot(Ss, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Trembling Aspen 
At<-subset(tree_dat, Species=="POPUTRE")
ggplot(At, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Paper birch 
Ep<-subset(tree_dat, Species=="BETUPAP")
ggplot(Ep, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)

#Ponderosa Pine
Py<-subset(tree_dat, Species=="PINUPON") 
ggplot(Py, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC) 

#White Spruce 
Sw<-subset(tree_dat, Species=="PICEGLA")  
ggplot(Sw, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC) 

#Western Larch 
Lw<-subset(tree_dat, Species=="LARIOCC") 
ggplot(Lw, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC) 

#Black Spruce 
Sb<-subset(tree_dat, Species=="PICEMAR") 
ggplot(Sb, aes(y=TotalA, x=year, color=MoistureRegime_clean))+
  geom_point()+
  facet_wrap(~GIS_BGC)



