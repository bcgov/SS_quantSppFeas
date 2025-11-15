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

#THIS SCRIPT DOES THE FOLLOWING: 
# pulls in cleaned tree data 
# cleans naming, crosswalks all plot data to newest BEC using classification list 
# For plots not indicated in classification list, used a combination of subzone and Mesoslope position, slope /aspect, Soil nutrient regime, 
# soil moisture regime and plot descriptions (depending on availability) to assign site series to plot data. 
# saves output file (tree_data_cleaned_updated.Rdata)

#libraries
library(tidyverse)


#crosswalk tree (plot) data---- 
#load BC tree data 
load(file="data/tree_data_cleaned.Rdata") 

#subset cols of interest 
tree_dat_sub<-dplyr::select(tree_dat, PlotNumber, Species, spp, TotalA, TotalB, 
                            UserSiteUnit, BECSiteUnit, 
                            GIS_BGC, 
                            SubZone,SiteSeries, #MapUnit, 
                            SitePlotQuality,NutrientRegime_clean,MoistureRegime_clean, 
                            Elevation, SlopeGradient, Aspect, MesoSlopePosition, Latitude, Longitude,
                            SuccessionalStatus, StructuralStage)   
str(tree_dat_sub)


#create matching columns to feas tables- from ss_cleaned.csv
#look at data across species/sites 
#match with updated site series info from WHM- BEC v12, v13 x plot numbers list
#ss_cleaned_old<-read.csv("data/All_BGC12DEC2024_SU.csv") 
ss_cleaned<-read.csv("data/All_BGC13_May2025_SU.csv") 

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
                       if_else(ss_nospace_final=="ESSFxv2/NA" & !is.na(UserSiteUnit), UserSiteUnit, ss_nospace_final))
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

tree_dat_sub<-subset(tree_dat_sub, !is.na(ss_nospace_final)& !grepl('NA', ss_nospace_final))#remove NA site series 

#save cleaned & updated tree data----
save(tree_dat_sub, file="data/tree_data_cleaned_updated.Rdata")

