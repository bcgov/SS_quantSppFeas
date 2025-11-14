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

#This script does the following:
#pulls in a csv of BEC plot data with no corresponding environmental suitability rating for a given spp x site series (from feas_tables.R)
#cross references to edatopic tables for matching site series
#quality filters plot data 
#assigns preliminary suitability rating using rule-set based on plot abundance data
#saves as output for further expert review and incorporation into suitability ratings dataset 

#libraries
library(tidyverse)

#read in plot data
BEC_nosuit<-read.csv("data/BEC_missing_feas.csv")%>%select(-X, -BECSiteUnit, -MapUnit)
BEC_nosuit<-mutate(BEC_nosuit, TotalA= if_else(is.na(TotalA), 0, TotalA),,TotalB= if_else(is.na(TotalB), 0, TotalB))
BEC_nosuit$TotalAB<-BEC_nosuit$TotalA+ BEC_nosuit$TotalB

#read in edatopic table 
edat<-read.csv("data/Edatopic_v13_11.csv")
edat<-filter(edat, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|SBAP|SASbo|MGPmg|MGPdm|
                          BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', SS_NoSpace))

ss_all<-unique(edat$SS_NoSpace) 
#subset by edat table
BEC_nosuit<-subset(BEC_nosuit, ss_nospace_final %in% ss_all)
ss_BEC<-unique(BEC_nosuit$ss_nospace_final) 

#calculate averages by site series by spp
avgs<-group_by(BEC_nosuit, ss_nospace_final, Species, spp, newsuit)%>%
  summarise(mean_abund_ss=mean(TotalAB, na.rm = T), sd_abund_ss=sd(TotalAB, na.rm = T), nplots_ss=n())%>%
  mutate(sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))

#quality filter 
avgsx<-subset(avgs, nplots_ss>2 & mean_abund_ss > sd_abund_ss) 

#assign suit ratings from plot averages
avgsx<-mutate(avgsx, newsuit= case_when(mean_abund_ss<1 ~ 4,
                                      mean_abund_ss>=1 & mean_abund_ss<10 ~ 3,
                                      mean_abund_ss>=10 & mean_abund_ss<25 ~ 2,
                                      mean_abund_ss>=24 ~ 1,
                                      TRUE~NA))
avgsy<-select(avgsx, ss_nospace_final, Species, spp, newsuit)
avgsy$bgc<-NULL

#join back with plot dataset
BEC_nosuit$newsuit<-NULL
BEC_nosuit<-left_join(BEC_nosuit, avgsy)

#select only infilled
BEC_nosuitx<-subset(BEC_nosuit, !is.na(newsuit))
rem<-anti_join(BEC_nosuit, BEC_nosuitx) #did not meet quality filter 

BEC_nosuitx$newsuit_ord<-as.ordered(BEC_nosuitx$newsuit)
BEC_nosuitx$ss_nospace<-NULL
BEC_nosuitx<-rename(BEC_nosuitx, ss_nospace=ss_nospace_final)
BEC_nosuitx$mod<-"imputed" 
BEC_nosuitx$bgc<-NULL
BEC_nosuitx<-separate(BEC_nosuitx, col =  ss_nospace, into = 'bgc',sep =  "/", remove = F)


BEC_nosuitx<-select(BEC_nosuitx,bgc,ss_nospace,suitability,spp,newsuit, 
                   mod,outrange, PlotNumber, Species,TotalA, TotalB,
                   NutrientRegime_clean, MoistureRegime_clean, Elevation,SlopeGradient,Aspect, MesoSlopePosition, 
                   Latitude, Longitude,SuccessionalStatus, StructuralStage,SiteUnit, newsuit_ord, TotalAB) 

write.csv(BEC_nosuitx, "data/imputed_suit_ratings.csv") 
#imp<-read.csv("data/imputed_suit_ratings.csv") 
