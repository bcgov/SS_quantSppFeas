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
#pulls in most recent feasibility tables 
#cleans naming, crosswalks to newest BEC 
#merges feas values with BC BEC plot data (tree abundance) & saves output (feasibility_abundance_data.Rdata)
#creates list of data for which we have BEC plot abundances but no feas ratings (BEC_missing_feas.csv)
#calculates and plots mean abundances by feas rating for top spp
#pulls in plot level climate data (from climr_getdata_plots.R)
#runs PCA on climate vars and merges back with feas & abundance values (sourced-climPCAs.R)
#saves rds file for further quantitative validation/modeling (feas_abund_clim_data.Rdata)

#ALL references to feasibility=> Environmental Suitability 

#libraries
library(tidyverse)

#load feasibility data---- 
feas.dat<-read.csv("data/Suitability_v13_19.csv") #v13_19 updated with Ecologist review as of May 15, 2025
feas.dat<-mutate(feas.dat, newsuit=if_else(is.na(newsuit), suitability, newsuit))#if not updated, use previous rating 

#clean and crosswalk tree (plot) data---- 
#take out the US and alberta stuff because it won't match plot data
feas.dat.sub<-filter(feas.dat, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))
#deal with duplicates for coastal and interior spp 
feas.dat.sub<-select(feas.dat.sub, -sppsplit) %>%distinct(.) #!="Plc"& sppsplit!="Pli"& sppsplit!="Fdc" & sppsplit!="Fdi"& sppsplit!="Pyc"& sppsplit!="Pyi")

#load BC tree data 
#load(file="data/tree_data_cleaned_wzeros.Rdata")  
load(file="data/tree_data_cleaned.Rdata") #only use actual plot data, not imputed zeros 

#subset cols of interest 
tree_dat_sub<-dplyr::select(tree_dat, PlotNumber, Species,TotalA, TotalB, 
                            UserSiteUnit, BECSiteUnit, GIS_BGC, 
                            SubZone,SiteSeries,        
                            MapUnit, SitePlotQuality,NutrientRegime_clean,MoistureRegime_clean, 
                            Elevation, SlopeGradient, Aspect, MesoSlopePosition, Latitude, Longitude,
                            SuccessionalStatus, StructuralStage)   
str(feas.dat.sub)
str(tree_dat_sub)

#add spp codes to plot data 
sort(unique(feas.dat.sub$spp))
sort(unique(tree_dat$Species))

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

nas<-subset(tree_dat_sub, is.na(spp))
unique(nas$Species)
rm(tree_dat)
tree_dat_sub<-subset(tree_dat_sub, !is.na(spp))
rm(nas)

#create matching columns to feas tables- from ss_cleaned.csv
#look at data across species/sites 
#match with updated site series info from WHM- BEC v12, v13 x plot numbers list
#ss_cleaned<-read.csv("data/All_BGC12DEC2024_SU.csv") 
ss_cleaned<-read.csv("data/All_BGC13_May2025_SU.csv") 
ss_cleaned$bgc<-NULL
#ss_cleaned<-rbind(ss_cleaned, ss_cleaned2)%>%distinct(.)
ss_cleaned<-mutate(ss_cleaned, SiteUnit= gsub("ESSFdv ", "ESSFdvw", SiteUnit))#dv not rated 
tree_dat_sub<-left_join(tree_dat_sub, ss_cleaned) 

tree_dat_sub<-mutate(tree_dat_sub,  ss_nospace= gsub(" ", "", SiteUnit))
tree_dat_sub<-separate(tree_dat_sub, col =  SiteUnit, into = 'bgc',sep =  "/", remove = F)%>%mutate(bgc=gsub(" ", "", bgc))


#deal with other data that not on ss_cleaned (~2/3 dataset)
tree_dat_sub<-mutate(tree_dat_sub, GIS_BGC= case_when(
  grepl('BWBS', GIS_BGC) & SubZone=='dk1'| SubZone=='dk 1'~'BWBSdk1',
  grepl('BWBS', GIS_BGC) & SubZone=='dk2'| SubZone=='dk 2'~'BWBSdk2',
  grepl('BWBS', GIS_BGC) & SubZone=='mw1'| SubZone=='mw 1'~'BWBSmw1',
  grepl('BWBS', GIS_BGC) & SubZone=='mw2'| SubZone=='mw 2'~'BWBSmw2',
  grepl('BWBS', GIS_BGC) & SubZone=='wk1'| SubZone=='wk 1'~'BWBSwk1',
  grepl('BWBS', GIS_BGC) & SubZone=='wk2'| SubZone=='wk 2'~'BWBSwk2',
  grepl('BWBS', GIS_BGC) & SubZone=='wk3'| SubZone=='wk 3'~'BWBSwk3',
  TRUE ~ GIS_BGC))

tree_dat_sub<-mutate(tree_dat_sub, SiteSeries = str_remove_all(SiteSeries, "ish|new|\\$|\\*|\\?|/|\\-|\\(|\\)|\\\\|YS|MS|k|z|p|u|FH|SH"))
tree_dat_sub<-mutate(tree_dat_sub, SiteSeries2=case_when(SiteSeries=="1"~"01",
                                                         SiteSeries=="2"~"02",
                                                         SiteSeries=="3"~"03",
                                                         SiteSeries=="4"~"04",
                                                         SiteSeries=="5"~"05",
                                                         SiteSeries=="6"~"06",
                                                         SiteSeries=="7"~"07",
                                                         SiteSeries=="8"~"08",
                                                         SiteSeries=="9"~"09",TRUE~NA))%>%
  mutate(SiteSeries=if_else(!is.na(SiteSeries2), SiteSeries2, SiteSeries))%>%select(-SiteSeries2)


tree_dat_sub<-unite(tree_dat_sub, col = "ss_nospace2", c("GIS_BGC", "SiteSeries"), sep = "/", remove=F)%>%
  relocate(ss_nospace2, .after = ss_nospace)

#remove symbols
tree_dat_sub<-mutate(tree_dat_sub,  ss_nospace=gsub("[$-]", "", ss_nospace))
tree_dat_sub<-mutate(tree_dat_sub,  ss_nospace2=gsub("[$-]", "", ss_nospace2))

#updated ss from LMH77 crosswalk 
cw<-read.csv("data/LMH77_BGC_su_Crosswalk_v3.csv") 
cw<-unite(cw, col = "ss_nospace2", c("OldBGC","OldSiteUnitLMH28_26"), sep = "/")
cw<-unite(cw, col = "ss_nospace_new", c("NewBGC","NewSiteUnit"), sep = "/")
cw<-select(cw, ss_nospace2, ss_nospace_new)%>%distinct(.)

#HAIDA GWAII
cw1<-read.csv("data/HaidaGwaii_cw.csv") 
cw1<-mutate(cw1, old.subzone= if_else(grepl('CWH', new.subzone), "CWH", "MH"))
cw1<-unite(cw1, col = "ss_nospace2", c("old.subzone","old.site.unit"), sep = "", remove=F)
cw1<-unite(cw1, col = "ss_nospace_new", c("new.subzone","new.site.unit"), sep = "/", remove=F)
cw1<-select(cw1, ss_nospace2, ss_nospace_new)%>%distinct(.)

cw<-rbind(cw, cw1)#combine HG & New coast b/c vh2 units overlap 

#BWBS update LMH65
#https://www2.gov.bc.ca/assets/gov/farming-natural-resources-and-industry/forestry/tree-species-selection/es_lists_bwbszone.xlsx
cw2<-read.csv('data/bwbs_crosswalk.csv')

#lmh75&76 crosswalk
cw3<-read.csv("data/LMH75_76_Crosswalk.csv") 

#ICH/ESSF-lmh 70&71 crosswalk 
cw4<-read.csv("data/LMH70_71cw.csv") 

#combine all crosswalk tables 
cw<-rbind(cw, cw2, cw3,cw4)

#join with BEC plot data 
tree_dat_sub<-left_join(tree_dat_sub, cw, relationship = "many-to-many")

#select appropriate column 
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=case_when(
  !is.na(ss_nospace)~ ss_nospace, #if on ss_cleaned, use
  is.na(ss_nospace)& !is.na(ss_nospace_new)~ ss_nospace_new, #if not on ss_cleaned, use cw 
  TRUE~ss_nospace2)) #if not on ss_cleaned or Cw, use use 'GIS_BGC, Site Series' cols

#other individual units/cases
#BG units 
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "964159"] <- "BGxw2/115"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "964119"] <- "BGxh3/111"	
tree_dat_sub$ss_nospace_final[tree_dat_sub$BECSiteUnit == "BG  xh 3 /113(Fm05)"] <- "BGxh3/113"	

#BWBS
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "BWBSdk/110"& tree_dat_sub$ss_nospace2=="BWBSdk/111"] <- "BWBSdk/111"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "966754"] <- "BWBSdk/110"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "BWBSwk3/NA" & tree_dat_sub$UserSiteUnit=="BWBSwk 3 /101-mixed"]<-"BWBSwk3/101"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "BWBSwk3/NA" & tree_dat_sub$UserSiteUnit=="BWBSwk 3 /102.1"]<-"BWBSwk3/102"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "92MK24"|tree_dat_sub$PlotNumber == "92MK22"] <- "BWBSmk/104b"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "BWBSdk/110"& tree_dat_sub$ss_nospace2=="BWBSdk/111"] <- "BWBSdk/111"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "BWBSmw/NA" & tree_dat_sub$MoistureRegime_clean==2]<-"BWBSmw/102"

#CDF/CWH
tree_dat_sub$ss_nospace_final[tree_dat_sub$UserSiteUnit == "CDF mm   /103Quergar-Camas"] <- "CDFmm/103" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="CWHvh2/09/10"] <- "CWHvh2/112"
tree_dat_sub$ss_nospace_final[tree_dat_sub$UserSiteUnit=="CWH xm 2 /08"] <- "CWHdm2/113"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="CWHwh1/103.1"& tree_dat_sub$spp!="Cw"] <- "CWHwh1/103" #only Cw rated on the variants 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="CWHws1/0607"] <- "CWHws1/113" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="CWHws1/0104"] <- "CWHws1/101" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_new=="CWHvh3/116"] <- "CWHvh3/116"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="CWHvh2/0910"] <- "CWHvh3/113"

#ESSF
tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="ESSFmm2/NA")%>%
  mutate(ss_nospace_new=case_when( edatope=="2C"~ "ESSFmm2/04", 
                                   edatope=="3A"~ "ESSFmm2/03", 
                                   edatope=="3B"~ "ESSFmm2/01", 
                                   edatope=="4C"~ "ESSFmm2/01", 
                                   edatope=="4B"~ "ESSFmm2/01", 
                                   edatope=="5C"~ "ESSFmm2/05", 
                                   edatope=="6B"~ "ESSFmm2/06", 
                                  TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="ESSFmm2/NA")%>%rbind(., tree_dat_subx)#bring back into main df 

tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="ESSFwv/NA"|ss_nospace_final=="ESSFwv/105"|ss_nospace_final=="ESSFwv/101")%>%
  mutate(ss_nospace_new=case_when( edatope=="1B"~ "ESSFwv/02", 
                                   edatope=="2B"| edatope=="3C"~ "ESSFwv/03", 
                                   edatope=="3B"~ "ESSFwv/04", 
                                   edatope=="4B"|edatope=="4C"~ "ESSFwv/01", 
                                   edatope=="4D"~ "ESSFwv/05", 
                                   edatope=="5C"|edatope=="5D"~ "ESSFwv/07", 
                                   edatope=="6D"~ "ESSFwv/06", 
                                   edatope=="7D"~ "ESSFwv/09", 
                                   TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="ESSFwv/NA")%>%rbind(., tree_dat_subx)#bring back into main df 

tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="ESSFvc/NA")%>%
  mutate(ss_nospace_new=case_when( edatope=="0B"|edatope=="0A" |edatope=="1B"|edatope=="1A"~ "ESSFvc/02", 
                                   edatope=="2C"|edatope=="2B"~ "ESSFvc/03", 
                                   edatope=="3C"|edatope=="3B"|edatope=="3D"| edatope=="4C"~ "ESSFvc/01", 
                                   edatope=="5C"| edatope=="5D"|edatope=="6D"~ "ESSFvc/04", 
                                   edatope=="6C"|edatope=="6B" ~ "ESSFvc/05", 
                                   edatope=="7D"|edatope=="7C"~ "ESSFvc/06", 
                                   TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="ESSFvc/NA")%>%rbind(., tree_dat_subx)#bring back into main df 

tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="ESSFxvw/NA")%>%
  mutate(ss_nospace_new=case_when( edatope=="2C"& MesoSlopePosition=="UP" ~ "ESSFxvw/02", 
                                   edatope=="2A"~ "ESSFxvw/03", 
                                   edatope=="2B"~ "ESSFxvw/04", 
                                   edatope=="3B"~ "ESSFxvw/06", 
                                   edatope=="4B"~ "ESSFxvw/01", 
                                   edatope=="6D"~ "ESSFxvw/09", 
                                   TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="ESSFxvw/NA")%>%rbind(., tree_dat_subx)#bring back into main df 

tree_dat_sub$ss_nospace_final[tree_dat_sub$UserSiteUnit == "MS mw/01"] <- "ESSFmw1/101"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "ESSFdvw/04"] <- "ESSFdvw/04"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace == "ESSFdvw/06"] <- "ESSFdvw/06"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "K001426"|tree_dat_sub$PlotNumber == "98-977"] <- "ESSFdvw/02"
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=
                       if_else(ss_nospace_final=="ESSFxv2/NA" & !is.na(BECSiteUnit), BECSiteUnit, ss_nospace_final))
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ESSFxv2/0304"] <- "ESSFxv2/03" 

tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(GIS_BGC=="ESSFxc3", ss_nospace2, ss_nospace_final))
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(grepl("ESSFdc3",ss_nospace_final), ss_nospace2, ss_nospace_final))
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ESSFxc2/103a"] <- "ESSFxc2/103"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ESSFxc2/NA"& tree_dat_sub$MoistureRegime_clean =='7' ] <- "ESSFxc2/112"
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(grepl("ESSFdv1|ESSFdv2",ss_nospace2), ss_nospace2, ss_nospace_final))

#ICH
tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="ICHmm/NA")%>%
  mutate(ss_nospace_new=case_when(edatope=="1A"| edatope=="1B"~ "ICHmm/02", 
                                  edatope=="2B"| edatope=="3B"~ "ICHmm/03", 
                                  edatope=="3C"|edatope=="4C"| edatope=="4B"~ "ICHmm/01", 
                                  edatope=="5C"|edatope=="5B"| edatope=="5D"~ "ICHmm/04",
                                  edatope=="6B"|edatope=="6C"| edatope=="6D"~ "ICHmm/05", 
                                  edatope=="6E"~ "ICHmm/06",
                                  edatope=="7B"~ "ICHmm/07",
                                  edatope=="7C"~ "ICHmm/08",
                                  TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="ICHmm/NA")%>%rbind(., tree_dat_subx)%>%select(-bgc)#bring back into main df 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="ICHdk/06"] <- "ICHdk/06" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="ICHdk/FF"] <- "ICHdk/04" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHdm/NA"& tree_dat_sub$MoistureRegime_clean =='1' ] <- "ICHdm/102"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHdw1/113(Ws10)"] <- "ICHdw1/113"
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(grepl("ICHdw3",ss_nospace_final), ss_nospace2, ss_nospace_final))
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "K12-089"] <- "ICHdw4/112"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "N770317"] <- "ICHmw2/111"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "U812023"] <- "ICHmw2/113"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "VRI682"|tree_dat_sub$PlotNumber == "VRI135"] <- "ICHmw4/101"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "VRI677"] <- "ICHdw4/101"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="ICHmw1/06"] <- "ICHmw1/06" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="ICHmw1/07"] <- "ICHmw1/07" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHmw1/102"] <- "ICHmw1/02" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHvk1/NA"& tree_dat_sub$MoistureRegime_clean =='0' ] <- "ICHvk1/02"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHvk1/NA"& tree_dat_sub$MoistureRegime_clean =='7' ] <- "ICHvk1/06"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="ICHwk1/NA"& tree_dat_sub$MoistureRegime_clean =='7' ] <- "ICHvk1/08"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace2=="ICHxwa/111"]<- "ICHxwa/111"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "694-C"] <- "ICHxwa/110"
tree_dat_sub$ss_nospace_final[tree_dat_sub$nospace_final=="ICHmw3/NA"& tree_dat_sub$MoistureRegime_clean =='7' ]<-"ICHmw3/09"


#IDF
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(grepl("IDFdc|IDFxc",ss_nospace2), ss_nospace2, ss_nospace_final))%>%
  mutate(ss_nospace_final=case_when(ss_nospace_final=="IDFxc/05a"~"IDFxc/05", ss_nospace_final=="IDFxc/05b"~"IDFxc/05", 
                                    ss_nospace_final=="IDFxc/08b"~"IDFxc/08", ss_nospace_final=="IDFxc/03b"~"IDFxc/03", TRUE~ss_nospace_final))
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="IDFmw2/101"& tree_dat_sub$ss_nospace2=="IDFmw2/05"] <- "IDFmw2/05"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_new=="IDFdm1/111"] <- "IDFdm1/111"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="PyGIF70"|tree_dat_sub$PlotNumber=="0109191"] <- "IDFdm1/110"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber == "VRI141"] <- "IDFdk2/101"


#MS
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final=if_else(grepl("MSdc1|MSxk3|MSdc3|MSdm3",ss_nospace2), ss_nospace2, ss_nospace_final))
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="MSxk3/0809"] <- "MSxk3/08" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="MSdc2/MP"] <- "MSdc2/08" 
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="91MK084"] <- "MSdv/03"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="91MK075"] <- "MSdv/07"
#SBPS

tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBPSdc/NA"& tree_dat_sub$MoistureRegime_clean =='6']<-"SBPSdc/08"

#SBS
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBSmm/102"] <- "SBSmm/02" 
tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="SBSdh2/NA")%>%
  mutate(ss_nospace_new=case_when(edatope=="1B"~"SBSdh2/03",
                                   edatope=="3A"~"SBSdh2/05",
                                  edatope=="5B"|edatope=="5C"~"SBSdh2/06",
                                  edatope=="2B"|edatope=="2C"~ "SBSdh2/04",
                                  edatope=="3B"| edatope=="3C"|edatope=="4C"| edatope=="4B"~ "SBSdh2/01", 
                                  edatope=="5D"|edatope=="6C"|edatope=="6D"|PlotNumber=="12-2317"~ "SBSdh2/07", 
                                  TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="SBSdh2/NA")%>%rbind(., tree_dat_subx)#bring back into main df

tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="8529210"] <- "SBSdw3/10"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="8102489"]<-"SBSmc1/03"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="8205791"]<-"SBSmc1/05"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="8529228"]<-"SBSmc2/08"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBSmc2/01A"]<-"SBSmc2/01a"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBSmc2/01C"]<-"SBSmc2/01ac"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBSmc2/01B"]<-"SBSmc2/01b"
tree_dat_sub$ss_nospace_final[tree_dat_sub$ss_nospace_final=="SBSmc2/0906"]<-"SBSmc2/09"
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="8102474"]<-"SBSmh/08"


tree_dat_subx<-unite(tree_dat_sub, edatope, MoistureRegime_clean, NutrientRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="SBSwk3a/NA")%>%
  mutate(ss_nospace_new=case_when(edatope=="3B"~"SBSwk3a/05",
                                  edatope=="3C"~"SBSwk3a/04",
                                  edatope=="4C"~"SBSwk3a/01",
                                  edatope=="5D"~"SBSwk3a/06",
                                  TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="SBSwk3a/NA")%>%rbind(., tree_dat_subx)#bring back into main df

tree_dat_subx<-unite(tree_dat_sub, edatope, NutrientRegime_clean, MoistureRegime_clean, sep="", remove=F)%>% 
  filter(ss_nospace_final=="SBSwk2/NA"|ss_nospace2=="SBSwk2/NA")%>%
  mutate(ss_nospace_new=case_when(edatope=="A6"~"SBSwk2/04", edatope=="B2"~"SBSwk2/02", 
                                  edatope=="B3"|edatope=="C3"~"SBSwk2/03",
                                  edatope=="C5"|edatope=="D5"~"SBSwk2/05",
                                  edatope=="C6"|edatope=="D6"~"SBSwk2/06",
                                  edatope=="C4"|edatope=="D3"~"SBSwk2/01",
                                  TRUE~NA))%>%
  mutate(ss_nospace_final=(if_else(!is.na(ss_nospace_new), ss_nospace_new,ss_nospace)))%>%select(-edatope)
tree_dat_sub<- filter(tree_dat_sub, ss_nospace_final!="SBSwk2/NA")%>%rbind(., tree_dat_subx)#bring back into main df


#ESSF
tree_dat_sub$ss_nospace_final[tree_dat_sub$PlotNumber=="O4CH036"] <- "ESSFmk/07"

#filter out poor quality 
tree_dat_sub<-subset(tree_dat_sub, !grepl("poor|Poor|POOR", SitePlotQuality))
tree_dat_sub<-subset(tree_dat_sub, !grepl("omit|Omit", UserSiteUnit))

#remove any stray spaces
tree_dat_sub<-mutate(tree_dat_sub, ss_nospace_final= gsub(" ", "", ss_nospace_final))

#correct un/unp/uns
tree_dat_subuns<-subset(tree_dat_sub, grepl('un|unp|uns|cdp|dkp|dvp|mcp|mkp|msp|mvp|mwp|mmp|vcp|wcp|wmp|whp|wvp|xcp|xvp', ss_nospace_final)) #pull out all site series 
tree_dat_sub<-anti_join(tree_dat_sub, tree_dat_subuns)#remove from main dataframe 
#set all to 00
tree_dat_subuns$add<-"00"
tree_dat_subuns<-separate(tree_dat_subuns,col = ss_nospace_final, into =  'ss_nospace_final', sep = '/')%>%
  unite(col = 'ss_nospace_final', c('ss_nospace_final', 'add'), sep = '/')
#update for x/m/h site series 
tree_dat_subuns2<-filter(tree_dat_subuns, grepl('BWBSvk|ESSFun|MHmmp|MHmsp|MHun|MHvhp|SBSun|SWBun|SWBvk', ss_nospace_final))%>%
  subset(ss_nospace_final!="ESSFunp/00" & ss_nospace_final!="MHunp/00"& ss_nospace_final!="SWBuns") #use 00
tree_dat_subuns<-anti_join(tree_dat_subuns, tree_dat_subuns2)
tree_dat_subuns2<-mutate(tree_dat_subuns2, add= case_when(MoistureRegime_clean=="0"|MoistureRegime_clean=="1"|MoistureRegime_clean=="2" ~"x", 
                        MoistureRegime_clean=="3"|MoistureRegime_clean=="4" ~"m",
                        MoistureRegime_clean=="5"|MoistureRegime_clean=="6"|MoistureRegime_clean=="7" ~"h"))%>%
  separate(col = ss_nospace_final, into =  'ss_nospace_final', sep = '/')%>%
  unite(col = 'ss_nospace_final', c('ss_nospace_final', 'add'), sep = '/')

tree_dat_subuns<-rbind(tree_dat_subuns, tree_dat_subuns2)
tree_dat_subuns<-subset(tree_dat_subuns, Elevation>900) #make sure it's all actually alpine or parkland 

tree_dat_sub<-rbind(tree_dat_sub, tree_dat_subuns)#bring back into main df
rm(tree_dat_subuns, tree_dat_subuns2)

#join BEC data with suitability data----
feas.dat.sub<-rename(feas.dat.sub, ss_nospace_final=ss_nospace)
feas.dat.subx<-left_join(feas.dat.sub, tree_dat_sub,  relationship = "many-to-many")  
feas.dat.suby<-subset(feas.dat.subx,!is.na(TotalA)| !is.na(TotalB))

feas.dat.subx2<-left_join(tree_dat_sub,  feas.dat.sub, relationship = "many-to-many")  
feas.dat.subz<-subset(feas.dat.subx2,!is.na(newsuit))

#set feas as ord factor
feas.dat.suby$newsuit_ord<-ordered(feas.dat.suby$newsuit, levels = c(5, 4,3, 2, 1))
str(feas.dat.suby$newsuit_ord)
knitr::kable(group_by(feas.dat.suby, newsuit_ord)%>%summarise(counts=n())) 

#remove duplicate cols
feas.dat.sub<-select(feas.dat.suby,  -UserSiteUnit,        
                     -BECSiteUnit, -GIS_BGC,  -SubZone, -SiteSeries, -MapUnit, -ss_nospace2, -ss_nospace_new,  -ss_nospace, 
                     -SitePlotQuality)%>%rename(ss_nospace =ss_nospace_final)   
feas.dat.sub<-distinct(feas.dat.sub)

#combine A & B layer
feas.dat.sub<-mutate(feas.dat.sub, TotalA= replace_na(TotalA, 0), TotalB= replace_na(TotalB, 0))%>%
  mutate(TotalAB=TotalA + TotalB)

save(feas.dat.sub, file="data/feasibility_abundance_data.Rdata") #save 

#identify things that may not have merged----
#have suit rating but no supporting plot data 
Feas_missing_BEC<-anti_join(feas.dat.subx, feas.dat.suby) 
Feas_missing_BEC<-subset(Feas_missing_BEC, suitability<5) #if E5 shouldn't expect in plot data 
x<-subset(Feas_missing_BEC, ss_nospace_final %in% tree_dat_sub$ss_nospace_final)#site series has plot data but not for all species 
y<-subset(tree_dat_sub, ss_nospace_final %in% Feas_missing_BEC$ss_nospace_final)
z<-anti_join(Feas_missing_BEC, x) #rated site series has no plot data for any spp 
missing_plots<-select(z, ss_nospace_final)%>%distinct(.) 

#sites with BEC data but missing suit ratings
BEC_missing_feas<-anti_join(feas.dat.subx2, feas.dat.subz)
BEC_missing_feas$X<-NULL
write.csv(BEC_missing_feas, "data/BEC_missing_feas.csv")

#feas.dat.suby<-rbind(feas.dat.suby, test)
#plot by spp ----
#create BGC zone column
feas.dat.sub<-mutate(feas.dat.sub, zone= case_when(grepl('ICH', ss_nospace)~"ICH",grepl('ESSF', ss_nospace)~"ESSF", grepl('MS', ss_nospace)~"MS",
                                grepl('SBPS', ss_nospace)~"SBPS", grepl('BAFA', ss_nospace)~"BAFA", grepl('CWH', ss_nospace)~"CWH", grepl('IDF', ss_nospace)~"IDF",
                                grepl('BG', ss_nospace)~"BG", grepl('ESSF', ss_nospace)~"ESSF", grepl('CDF', ss_nospace)~"CDF", grepl('SBS', ss_nospace)~"SBS", 
                                grepl('MH', ss_nospace)~"MH", grepl('CMA', ss_nospace)~"CMA", grepl('PP', ss_nospace)~"PP", grepl('BWBS', ss_nospace)~"BWBS", TRUE~ NA))

#calculate average abundances by feas scores
avgs<-group_by(feas.dat.sub, zone, bgc, ss_nospace, Species, spp, newsuit_ord)%>%
  summarise(mean_abund_ss=mean(TotalAB, na.rm = T), sd_abund_ss=sd(TotalAB, na.rm = T), nplots_ss=n())
  #%>%ungroup(.)%>%
  #group_by(bgc)#%>% mutate(mean_abund_subzone=mean(mean_abund_ss, na.rm = T))%>% ungroup(.)%>%
  #group_by(zone)%>% mutate(mean_abund_zone=mean(mean_abund_subzone, na.rm = T))
sort(unique(avgs$spp))
#"Ac" "At" "Ba" "Bg" "Bl" "Cw" "Dr" "Ep" "Fd" "Hm" "Hw" "La" "Lt" "Lw" "Mb" "Mv" "Pa" "Pl" "Pw" "Py" "Qg" "Ra" "Sb"
#"Ss" "Sx" "Tw" "Yc"
avgs<-mutate(avgs, sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))

#calculate potential cutoffs requiring additional review
#avgs<-group_by(avgs, newsuit_ord, spp, zone)#%>%mutate(higher=quantile(mean_abund_ss, 0.75),
                                            #          lower=quantile(data, 0.25))%>%
                                            #          mutate(max=lower-1.5*IQR(mean_abund_ss),
                                            #                 min(higher+1.5*IQR(mean_abund_ss)))

#spp plots 
#Ac- 
ggplot(subset(avgs, spp=="Ac"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  #geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("At")
#At- 
ggplot(subset(avgs, spp=="At"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  #geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("At")
#Ba- 
ggplot(subset(avgs, spp=="Ba"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Ba")
#Bg-
ggplot(subset(avgs, spp=="Bg"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Bg")
#Bl- 
ggplot(subset(avgs, spp=="Bl"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Bl")
#Cw-  
ggplot(subset(avgs, spp=="Cw"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Cw")
#Dr- 
ggplot(subset(avgs, spp=="Dr"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ggtitle("Dr")
#Ep-
ggplot(subset(avgs, spp=="Ep"), aes(x = newsuit_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ggtitle("Ep")
#Fd - 
ggplot(subset(avgs, spp=="Fd"), aes(x = newsuit_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Fd")+
  xlab("Environmental Suitability") + ylab("Relative abundance (% cover)") + xlab("Environmental Suitability")
#Hm- 
ggplot(subset(avgs, spp=="Hm"), aes(x = newsuit_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hm")
#Hw- 
ggplot(subset(avgs, spp=="Hw"), aes(x = newsuit_ord, y = mean_abund_ss,, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hw")
#Lw-
ggplot(subset(avgs, spp=="Lw"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Lw")
#Mb-
ggplot(subset(avgs, spp=="Mb"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Mb")
#Pa- 
ggplot(subset(avgs, spp=="Pa"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pa")
#Pl- 
ggplot(subset(avgs, spp=="Pl"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pl")
#Pw-  
ggplot(subset(avgs, spp=="Pw"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pw")
#Py- 
ggplot(subset(avgs, spp=="Py"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Py")
#Ra- 
ggplot(subset(avgs, spp=="Ra"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ra")
#Sb- 
ggplot(subset(avgs, spp=="Sb"), aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Sb")
#Sx-
ggplot(subset(avgs, spp=="Sx"), aes(x = newsuit_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Sx")
#Ss-
ggplot(subset(avgs, spp=="Ss"), aes(x = newsuit_ord, y = mean_abund_ss,alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ss")
#Yc - 
ggplot(subset(avgs, spp=="Yc"), aes(x = newsuit_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Yc")


#pull in climate data ----
##select climate variables from BGC model
climrVars = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
              "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
              "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
              "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", "PPT_JAS", "CMD.total")
names(feas.dat) 

#from climr_getdata_plots.R
load("data/clim_dat.plots.Rdata")
clim_dat<-rename(clim_dat, PlotNumber=id)%>% select(c("PlotNumber", "PERIOD"), climrVars)

feas.dat.clim<-left_join(feas.dat.sub, clim_dat)

save(feas.dat.clim, file="data/feas_abund_clim_data.Rdata")

#calculate PC axes for climate params----
source("scripts/Bayesian_analysis/climPCAs.R")


##MISC----

#other plot data-right now not at SS level so can't merge with feas tables #### 12/2024
#AB data 
#load US tree data 
#load( file= )
#names(feas.dat.sub)
#cols<-names(feas.dat.sub)
#all_US_dat<-select(all_US_dat)
#all_US_dat<-subset(all_US_dat, select = names(all_US_dat) %in% cols) 

#pulling in from Will's Ordinal forest model script 
feas.dat.clim <- readRDS("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_workingfiles/Feasibility_modelling/OrdinalForest_data.rds")
feas.dat.clim<-na.omit(feas.dat.clim)

feas.dat<-mutate(feas.dat, CMD.total=CMD.def +CMD)
varsl<-c(c("zone", "bgc", "ss_nospace","spp", "newfeas",  "fid", "WNA_DEM_4326_clipped" ,"xcoord", "ycoord") , climrVars)
feas.dat.sub<-select(feas.dat, varsl)

names(feas.dat)
names(feas.dat.clim)
feas.dat.clim<-rename(feas.dat.clim, bgc=BGC, ss_nospace=SS_NoSpace)

#remove previous feasibility cols 
feas.dat.clim<-select(feas.dat.clim, -newfeas, -ESuit)
feas.dat$feasible<-NULL

#join feasibility (enviro suitability) w/ climate data 
feas.dat<-left_join(feas.dat, feas.dat.clim)

#subset to only top 16 spp
spp_tab0<-tree_dat%>%  group_by(Species)%>%  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
feas_tab<-subset(feas_tab, Species %in% spp_keep$Species)
rm(spp_tab0)


