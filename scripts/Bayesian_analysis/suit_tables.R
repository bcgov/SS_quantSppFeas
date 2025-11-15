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
#subsets for BC feas ratings 
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
#feas.dat<-read.csv("data/Suitability_v13_19.csv") #v13_19 updated with Ecologist review as of May 15, 2025
feas.dat<-read.csv("data/Suitability_v13_22.csv") #v13_22 updated with Craig Delong review & inputed ratings Oct 1, 2025
feas.dat<-subset(feas.dat, spp!='X')#remove any with no species defined 

feas.dat<-mutate(feas.dat, newsuit=if_else(is.na(newsuit), suitability, newsuit))#if not updated, use previous rating 


#take out the US and alberta stuff because it won't match plot data
feas.dat.sub<-filter(feas.dat, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))

#load BEC plot data
load(file="data/tree_data_cleaned_updated.Rdata") 
#join BEC data with suitability data----
feas.dat.sub<-rename(feas.dat.sub, ss_nospace_final=ss_nospace)
feas.dat.subx<-left_join(feas.dat.sub, tree_dat_sub,by = c('ss_nospace_final', 'spp'),relationship = "many-to-many")  
feas.dat.suby<-subset(feas.dat.subx,!is.na(TotalA)| !is.na(TotalB))

feas.dat.subx2<-left_join(tree_dat_sub,  feas.dat.sub, relationship = "many-to-many")  
feas.dat.subz<-subset(feas.dat.subx2,!is.na(newsuit))

#set feas as ord factor
feas.dat.suby$newsuit_ord<-ordered(feas.dat.suby$newsuit, levels = c(5, 4,3, 2, 1))
str(feas.dat.suby$newsuit_ord)
knitr::kable(group_by(feas.dat.suby, newsuit_ord)%>%summarise(counts=n())) 

#remove duplicate cols
feas.dat.sub<-select(feas.dat.suby,  -UserSiteUnit,        
                     #-BECSiteUnit,-MapUnit,
                     -GIS_BGC,  -SubZone, -SiteSeries,  -ss_nospace2, -ss_nospace_new,  -ss_nospace, 
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

#sites with BEC data but missing suit ratings- for inputting 
BEC_missing_feas<-anti_join(feas.dat.subx2, feas.dat.subz)
BEC_missing_feas$X<-NULL
#write.csv(BEC_missing_feas, "data/BEC_missing_feas.csv")

#plot by spp ----
#create BGC zone column
feas.dat.sub<-mutate(feas.dat.sub, zone= case_when(grepl('ICH', ss_nospace)~"ICH",grepl('ESSF', ss_nospace)~"ESSF", grepl('MS', ss_nospace)~"MS",
                                grepl('SBPS', ss_nospace)~"SBPS", grepl('BAFA', ss_nospace)~"BAFA", grepl('CWH', ss_nospace)~"CWH", grepl('IDF', ss_nospace)~"IDF",
                                grepl('BG', ss_nospace)~"BG", grepl('ESSF', ss_nospace)~"ESSF", grepl('CDF', ss_nospace)~"CDF", grepl('SBS', ss_nospace)~"SBS", 
                                grepl('MH', ss_nospace)~"MH", grepl('CMA', ss_nospace)~"CMA", grepl('PP', ss_nospace)~"PP", grepl('BWBS', ss_nospace)~"BWBS", TRUE~ NA))

#calculate average abundances by feas scores
avgs<-group_by(feas.dat.sub, zone, bgc, ss_nospace, Species, spp, newsuit_ord)%>%
  summarise(mean_abund_ss=mean(TotalAB, na.rm = T), sd_abund_ss=sd(TotalAB, na.rm = T), nplots_ss=n())
avgs<-mutate(avgs, sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))

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

