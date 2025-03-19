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


#read in and clean/collate US & AB plot data 
veg_dat<-read.csv("C:/Users/ccollins/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/Feasibility_modelling/USABEC_Veg.csv")
plot_dat<-read.csv("C:/Users/ccollins/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/Feasibility_modelling/USABEC_Env.csv")

all_US_dat<-left_join(veg_dat, plot_dat, by="PlotNumber")#join
# no site series calls in this data ???

names(all_US_dat)

# SNR, aSMR 
unique(all_US_dat$NutrientRegime)# not recorded 
unique(all_US_dat$MoistureRegime) #looks good 

#make another col so will join with BEC data
all_US_dat$MoistureRegime_clean<-all_US_dat$MoistureRegime

#all_US_dat<-rename(zone=Zone, )

#AB 
#pull in data 
plot_datAB<-read.csv("C:/Users/ccollins/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/Feasibility_modelling/AllAlbertaEcoData/Plots.csv")
veg_datAB<-read.csv("C:/Users/ccollins/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/Feasibility_modelling/AllAlbertaEcoData/VegetationSpecies.csv")
veg_dat2AB<-read.csv("C:/Users/ccollins/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/Feasibility_modelling/AllAlbertaEcoData/VegetationTotals.csv")
veg_datAB<-left_join(veg_datAB, veg_dat2AB)

#reorg 
names(BEC_data$veg) #want to make column structure consistent 
#fix spaces in naming 
veg_datAB <- janitor::clean_names(veg_datAB)
#remove cols with no data 
veg_datAB<-janitor::remove_empty(veg_datAB, which = "cols")
names(veg_datAB)

#AB data -rename cols that already exist
veg_datAB<-rename(veg_datAB, PlotNumber=field_plot_number, Species=species_code, Cover2= overstory_tree_percent_cover,Height2=overstory_tree_height,
                Cover3=understory_tree_percent_cover, Height3=understory_tree_height, Cover4=tall_shrub_percent_cover, Height4=tall_shrub_height, 
                Cover5=low_shrub_percent_cover, Height5=low_shrub_height, Cover10=epiphytes_percent_cover, Layer=stratum_code)
#create 
veg_datAB<-mutate(veg_datAB, TotalA = Cover2 + Cover3)%>%mutate(HeightA = Height2 + Height3)%>%mutate(TotalB = Cover4 + Cover5)%>%
  mutate(HeightB = Height4 + Height5)%>%mutate(Cover6 = low_forb_percent_cover + graminoid_percent_cover)%>%mutate(Height6 = low_forb_height+ graminoid_height)%>%
  mutate(Cover7 = lichen_percent_cover + moss_percent_cover)

veg_datAB$ID<-row.names(veg_datAB)

#AB only uses 2 tree classes not 3, so put overstory into 2 and understory into 3 
veg_datAB$Cover1<-NA
veg_datAB$Height1<-NA

#add additional required cols from BEC db - these are all blank or inconsistently reported in BEC db 
veg_datAB$Cover5a<-NA
veg_datAB$Cover5b<-NA
veg_datAB$Cover5c<-NA
veg_datAB$Height5a<-NA
veg_datAB$Height5b<-NA
veg_datAB$Height5c<-NA
veg_datAB$Cover8<-NA
veg_datAB$Cover9<-NA
veg_datAB$LL<-NA
veg_datAB$AF<-NA
veg_datAB$DC<-NA
veg_datAB$UT<-NA
veg_datAB$VI<-NA
veg_datAB$PV<-NA
veg_datAB$PG<-NA
veg_datAB$FFA<-NA
veg_datAB$Cultural1<-NA
veg_datAB$Cultural2<-NA
veg_datAB$Other1<-NA
veg_datAB$Other2<-NA
veg_datAB$Collected<-NA
veg_datAB$Flag<-NA

#select and reorder cols to keep 
veg_datAB<-select(veg_datAB, PlotNumber, Species, Layer, Cover1, Height1,  Cover2,  Height2, Cover3,  Height3, TotalA,   HeightA ,Cover4,    
Height4,    Cover5,     Height5,    Cover5a,    Height5a,   Cover5b,    Height5b,   Cover5c,    Height5c,   TotalB,     HeightB,    Cover6,    
Height6,    Cover7,     Cover8,     Cover9,     Cover10,    Collected,  Flag,       ID,         LL,         AF,         DC,         UT,        
 VI,         PV,         PG,         FFA,        Cultural1,  Cultural2,  Other1,     Other2)    
