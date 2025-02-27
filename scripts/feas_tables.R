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
feas.dat<-read.csv("data/Feasibility_v13_4.csv")
feas.dat<-subset(feas.dat, !is.na(newfeas))

#merge feas & tree data---- 

#take out the US and alberta stuff because it won't match plot data
feas.dat.sub<-filter(feas.dat, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|SBAP|SASbo', ss_nospace))

#deal with duplicates for coastal and interior spp 
feas.dat.sub<-select(feas.dat.sub, -sppsplit) %>%distinct(.) #!="Plc"& sppsplit!="Pli"& sppsplit!="Fdc" & sppsplit!="Fdi"& sppsplit!="Pyc"& sppsplit!="Pyi")

#load BC tree data 
#load(file="data/tree_data_cleaned_wzeros.Rdata")  
load(file="data/tree_data_cleaned.Rdata") #only use actual plot data, not imputed zeros 

#subset cols of interest 
tree_dat_sub<-dplyr::select(tree_dat, PlotNumber, Species,TotalA, SiteUnit, Site, bgc, NutrientRegime_clean,MoistureRegime_clean, 
                            Elevation, SlopeGradient, Aspect, MesoSlopePosition, Latitude, Longitude)
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
                                                      Species=="ALNURUB"~"Dr",
                                                      Species=="BETUPAP"~"Ep", 
                                                      Species=="CALLNOO"~"Yc", 
                                                      Species=="LARIOCC"~"Lw",  
                                                      Species=="PICEENE"~ "Sx", #not differentiated
                                                      Species=="PICEGLA"~"Sx", #not differentiated
                                                      Species=="PICEMAR"~"Sb",
                                                      Species=="PICESIT"~"Ss",#not differentiated
                                                      Species=="PICEXLU"~"Ss",#not differentiated
                                                      Species=="PINUALB"~"Pa",
                                                      Species=="PINUCON"~"Pl",#not differentiated
                                                      Species=="PINUPON"~"Py",
                                                      Species=="PINUMON"~"Pw",
                                                      Species=="POPUTRE" ~"At",
                                                      Species=="POPUBAL"~"Ac", #not differentiated
                                                      Species=="POPUTRI"~"Act", 
                                                      Species=="PSEUMEN"~"Fd", #not differentiated 
                                                      Species=="TAXUBRE" ~"Tw",
                                                      Species=="THUJPLI"~"Cw",
                                                      Species=="TSUGHET" ~"Hw",
                                                      Species=="TSUGMER"~"Hm")) 




#update site series naming in plot data
#ss <- read.csv("C:/Users/ccollins/OneDrive - Government of BC/CCISS/ccissv13_workingfiles/Feasibility_modelling/All_BGC12DEC2024_SU.csv") #need to update to BEC 13 list when avaialable 
#ss<-mutate(ss, bgc= gsub(" ", "", bgc))%>%separate(SiteUnit, into = "zone", sep = " ", remove=F)%>%
#       mutate(zone2= case_when(grepl("ESSF", zone)~"ESSF",
#                                grepl("BWBS", zone)~"BWBS",
#                               grepl("SBPS", zone)~"SBPS", TRUE~NA))%>%mutate(zone=if_else(!is.na(zone2), zone2, zone))%>%select(-zone2)
#tree_dat_sub<-rename(tree_dat_sub, SiteUnitold=SiteUnit)
#tree_dat_sub<-left_join(tree_dat_sub, ss)
#tree_dat_sub<-mutate(tree_dat_sub, ss_new= if_else(SiteUnitold==SiteUnit, 'N', 'Y'))

#create matching columns to feas tables
tree_dat_sub<-mutate(tree_dat_sub,  ss_nospace= gsub(" ", "", SiteUnit))
tree_dat_sub<-mutate(tree_dat_sub,  bgc= gsub(" ", "", bgc))

#join feas with BEC data
feas.dat.subx<-left_join(feas.dat.sub, tree_dat_sub,  relationship = "many-to-many")  
feas.dat.suby<-subset(feas.dat.subx,!is.na(TotalA))


feas.dat.subx2<-left_join( tree_dat_sub,  feas.dat.sub, relationship = "many-to-many")  
feas.dat.subz<-subset(feas.dat.subx2,!is.na(newfeas))

#create list with missing feas ratings for imputing
BEC_missing_feas<-anti_join(feas.dat.subx2, feas.dat.subz)
Feas_missing_BEC<-anti_join(feas.dat.subx, feas.dat.suby)

sort(unique(BEC_missing_feas$ss_nospace))
BEC_missing_feas_unique<-select(BEC_missing_feas, ss_nospace, spp) %>% distinct(.)
BEC_missing_feas_unique<-filter(BEC_missing_feas_unique, !grepl('support', ss_nospace)) 
unique(BEC_missing_feas_unique$ss_nospace)
write.csv(BEC_missing_feas_unique, "data/BEC_missing_feas.csv")

#set feas as ord factor
feas.dat.suby$newfeas_ord<-ordered(feas.dat.suby$newfeas, levels = c(5, 4,3, 2, 1))
str(feas.dat.suby$newfeas_ord)
knitr::kable(group_by(feas.dat.suby, newfeas_ord)%>%summarise(counts=n())) #30 E5s with non-zero abundances 

#save feas plus plot data
feas.dat.sub<-feas.dat.suby

save(feas.dat.sub, file="data/feasibility_abundance_data.Rdata")


#plot by spp with plot data 
#create BGC zone column
feas.dat.sub<-separate(feas.dat.sub, SiteUnit, into = "zone", sep = " ", remove=F)%>%
      mutate(zone2= case_when(grepl("ESSF", zone)~"ESSF",
                                  grepl("BWBS", zone)~"BWBS",
                                 grepl("SBPS", zone)~"SBPS", TRUE~NA))%>%mutate(zone=if_else(!is.na(zone2), zone2, zone))%>%select(-zone2)
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
#Ac- fine- rare- needs more bgcs?
ggplot(subset(avgs, spp=="Ac"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ac")
#At- too high 2s and 3s
ggplot(subset(avgs, spp=="At"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("At")
#Ba- looks good 
ggplot(subset(avgs, spp=="Ba"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Ba")
#Bg- CWH 2s high
ggplot(subset(avgs, spp=="Bg"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Bg")
#Bl- ESSF- some high 2s and 3s 
ggplot(subset(avgs, spp=="Bl"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ ggtitle("Bl")
#Cw- high 3s in IDF & ESSF
ggplot(subset(avgs, spp=="Cw"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Cw")
#Dr- looks ok - some high 2s in CWH
ggplot(subset(avgs, spp=="Dr"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ggtitle("Dr")
#Ep- high 2s and 3s 
ggplot(subset(avgs, spp=="Ep"), aes(x = newfeas_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')+ggtitle("Ep")
#Fd - looks good 
ggplot(subset(avgs, spp=="Fd"), aes(x = newfeas_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Fd")
#Hm- looks ok- some high 2s and 3s in CWH
ggplot(subset(avgs, spp=="Hm"), aes(x = newfeas_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hm")
#Hw- looks good  
ggplot(subset(avgs, spp=="Hw"), aes(x = newfeas_ord, y = mean_abund_ss,, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Hw")
#Lw-looks good, high 2s in ICH 
ggplot(subset(avgs, spp=="Lw"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Lw")
#Mb-looks good, high 2s in CWH 
ggplot(subset(avgs, spp=="Mb"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Mb")
#Pa- not good- 3s higher than 2
ggplot(subset(avgs, spp=="Pa"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pa")
#Pl- looks good, some high 2s in CWH and IDF
ggplot(subset(avgs, spp=="Pl"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pl")
#Pw- looks good, some high 2s in CWH 
ggplot(subset(avgs, spp=="Pw"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Pw")
#Py- looks good
ggplot(subset(avgs, spp=="Py"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Py")
#Ra- very rare 
ggplot(subset(avgs, spp=="Ra"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ra")
#Sb- BWBS- high 3s and high 2s in SBS
ggplot(subset(avgs, spp=="Sb"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Sb")
#Sx-looks good- high E2 high in ESSF, IDF
ggplot(subset(avgs, spp=="Sx"), aes(x = newfeas_ord, y = mean_abund_ss,  alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Sx")
#Ss- high E4s in MH
ggplot(subset(avgs, spp=="Ss"), aes(x = newfeas_ord, y = mean_abund_ss,alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none') +ggtitle("Ss")
#Yc - good, some high 3s in CWH
ggplot(subset(avgs, spp=="Yc"), aes(x = newfeas_ord, y = mean_abund_ss,  alpha=0.5))+
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


