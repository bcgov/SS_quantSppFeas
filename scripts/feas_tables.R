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
#pulls in most recent feasibility tables (v 13_3) 
#cleans naming
#merges feas values with BC BEC plot data (tree abundance)
#calculates and plots mean abundances by feas rating for top spp
#pulls in plot level climate data (from climr_getdata_plots.R)
#runs PCA on climate vars and merges back with feas & abundance values (sourced-climPCAs.R)
#creates list of data for which we have BEC plot abundances but no feas ratings (BEC_missing_feas)- need to cross check by site series 
#saves rds file for further quantitative validation/modeling (feas_abund_clim_data.Rdata)

#libraries
library(tidyverse)

#load feasibility data---- 
feas.dat<-read.csv("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_latest_tool_materials/Feasibility_v13_3.csv")
feas.dat<-subset(feas.dat, !is.na(newfeas))

#filter out small BGCs not included in training/projections (from BGC projections repo- create training set.R)- TBD 2/12/25- waiting on updated BGC projections from Colin 
#BGC_exclude<-c( "BWBScm"  ,"CMAun_OR",  "CMAwh",   "CWHxs",  "ESSFdcp", "ESSFdh1" ,"ESSFdh2"  ,"ESSFdh_WA","ESSFwh_MT","ESSFxcp", "ESSFxcw",  
#                "ESSFxh_WA", "ESSFxvw" ,  "ICHmc1a",  "ICHxwa"  ,  "IDFww"  ,   "IDFww1"   , "IDFww2"  ,  "IDFxx1"  ,  "MHun"  ,    "MHunp"  ,   "MHwh" ,    
#                "MHwhp"  ,   "MSdc2"  ,   "MSun")

#feas.dat.sub<-subset(feas.dat, !(BGC %in% BGC_exclude))


#merge feas & tree data---- 
#deal with duplicates for coastal and interior spp 
feas.dat.sub<-subset(feas.dat, spp!="Plc"& spp!="Pli"& spp!="Fdc" & spp!="Fdi"& spp!="Pyc"& spp!="Pyi")

#add full spp codes to feas df 
feas.dat.subx<-mutate(feas.dat.sub, Species= case_when(spp=="Ba"~"ABIEAMA",
                                                      spp=="Ra" ~"ARBUMEN",
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
                                                     # spp=="Plc"~"PINUCON", #coastal
                                                    #  spp=="Pli"~"PINUCON", #interior
                                                    #  spp=="Pyc"~"PINUPON",
                                                    #  spp=="Pyi"~"PINUPON",
                                                      spp=="Py"~"PINUPON",
                                                      spp=="Pw"~"PINUMON",
                                                      spp=="Acb"~"POPUBAL",
                                                      spp=="At"~"POPUTRE",
                                                      spp=="Act"~"POPUTRI",
                                                      spp=="Fd"~"PSEUMEN", 
                                                     # spp=="Fdi"~"PSEUMEN", 
                                                     # spp=="Fdc"~"PSEUMEN", 
                                                      spp=="Tw"~ "TAXUBRE",
                                                      spp=="Cw"~ "THUJPLI",
                                                      spp=="Hw"~ "TSUGHET",
                                                      spp=="Hm"~"TSUGMER")) 


#take out the US and alberta stuff because it won't match plot data
feas.dat.subx<-filter(feas.dat.subx, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abE|abC|SBAP|SASbo|PPxh|MSd|MSx', ss_nospace))

#load BC tree data 
load(file="data/tree_data_cleaned_wzeros.Rdata")  
#load(file="data/tree_data_cleaned.Rdata")
tree_dat<-tree_dat_wzeros
rm(tree_dat_wzeros)
gc()

#update site series naming plot data
ss <- read.csv("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_workingfiles/Feasibility_modelling/All_BGC12DEC2024_SU.csv")
ss<-mutate(ss, bgc= gsub(" ", "", bgc))%>%separate(SiteUnit, into = "zone", sep = " ", remove=F)%>%
       mutate(zone2= case_when(grepl("ESSF", zone)~"ESSF",
                                grepl("BWBS", zone)~"BWBS", TRUE~NA))%>%mutate(zone=if_else(!is.na(zone2), zone2, zone))%>%select(-zone2)
tree_dat<-rename(tree_dat, SiteUnitold=SiteUnit)
tree_dat<-left_join(tree_dat, ss)
tree_dat<-mutate(tree_dat, ss_new= if_else(SiteUnitold==SiteUnit, 'N', 'Y'))

#subset cols of interest 
tree_dat_sub<-dplyr::select(tree_dat, PlotNumber, Species, TotalA, SiteUnitold, ss_new,
                     SiteUnit, bgc, zone, NutrientRegime_clean,MoistureRegime_clean, Latitude, Longitude)

#create matching columns to feas tables
tree_dat_sub<-mutate(tree_dat_sub,  ss_nospace= gsub(" ", "", SiteUnit))

#join with BC data
feas.dat.suby<-left_join(feas.dat.subx, tree_dat_sub)  
feas.dat.suby<-subset(feas.dat.suby,!is.na(TotalA))

#create list with missing feas ratings for imputing
feas.dat.subz<-left_join(tree_dat_sub, feas.dat.subx)
BEC_missing_feas<-subset(feas.dat.subz, is.na(newfeas))

#set feas as ord factor
feas.dat.suby$newfeas_ord<-ordered(feas.dat.suby$newfeas, levels = c(5, 4,3, 2, 1))
str(feas.dat.suby$newfeas_ord)
knitr::kable(group_by(feas.dat.suby, newfeas_ord)%>%summarise(counts=n())) #4 is actually a 5 here (i.e. zero)

#save feas plus plot data
feas.dat.sub<-feas.dat.suby

save(feas.dat.sub, file="data/feasibility_abundance_data.Rdata")


#plot by spp with plot data 
#calculate average abundances by feas scores
avgs<-group_by(feas.dat.sub, zone, bgc, ss_nospace, Species, spp, newfeas_ord)%>%summarise(mean_abund_ss=mean(TotalA, na.rm = T), sd_abund_ss=sd(TotalA, na.rm = T))%>%ungroup(.)%>%
  group_by(bgc)%>% mutate(mean_abund_subzone=mean(mean_abund_ss, na.rm = T))%>% ungroup(.)%>%
  group_by(zone)%>% mutate(mean_abund_zone=mean(mean_abund_subzone, na.rm = T))
avgs<-na.omit(avgs)

#how many spp have >1 rating?
check<-group_by(feas.dat.sub, spp)%>%summarise(n_ratings= n_distinct(newfeas), abund=mean(TotalA, na.rm=T))
check<-subset(check, n_ratings>1 & !is.na(abund))
check$spp

#spp plots 
#At- fine, some high 2s and 3s
ggplot(subset(avgs, Species=="POPUTRE"), aes(x = newfeas_ord, y = log(mean_abund_ss+1), colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("At")
#Ba- good, rare
ggplot(subset(avgs, Species=="ABIEAMA"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ba")
#Bl- CWH, IDF, MH look off 
ggplot(subset(avgs, Species=="ABIELAS"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Bl")
#Cw- issues in IDF & MH
ggplot(subset(avgs, Species=="THUJPLI"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Cw")
#Ep- looks good, some high 2s and 3s 
ggplot(subset(avgs, Species=="BETUPAP"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ggtitle("Ep")
#Fd - looks good 
ggplot(subset(avgs, Species=="PSEUMEN"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Fd")
#Hm- fine, rare 
ggplot(subset(avgs, Species=="TSUGMER"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hm")
#Hw- looks good  
ggplot(subset(avgs, Species=="TSUGHET"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hw")
#Lw-looks good, rare 
ggplot(subset(avgs, Species=="LARIOCC"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Lw")
#Pl- good
ggplot(subset(avgs, Species=="PINUCON"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pl")
#Se- fine, rare, E2 high in SBS
ggplot(subset(avgs, Species=="PICEMAR"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Se")
#Ss- good
ggplot(subset(avgs, Species=="PICESIT"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ss")
#Yc - good, need more ratings in ESSF? 
ggplot(subset(avgs, Species=="CALLNOO"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Yc")

#look at E5s with non-zero abundances
check<-subset(avgs, newfeas_ord=="5" & mean_abund_ss>0) #very few 


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
source("scripts/climPCAs.R")

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

#save again
save(feas.dat.sub, file="data/feasibility_data.Rdata")

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


