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
#pulls in most recent suitability tables 
#subsets for BC suit ratings 
#merges suit ratings with BC BEC plot data (tree abundance) & saves output (feasibility_abundance_data.Rdata)
#creates list of data for which we have BEC plot abundances but no feas ratings (BEC_missing_feas.csv)
#calculates and plots mean abundances by feas rating for top spp
#quality filters so that all suitability ratings had a minimum of 2 plots for validation 
#creates a unique variable for the 15 edatopic spaces for modeling 
#saves dataset for modeling validation 
#OPTIONAL 
#pulls in plot level climate data (from climr_getdata_plots.R)
#runs PCA on climate vars and merges back with feas & abundance values (sourced-climPCAs.R)
#saves rds file for further quantitative validation/modeling (feas_abund_clim_data.Rdata)


#ALL references to feasibility=> Environmental Suitability 

#libraries
library(tidyverse)

#load feasibility data---- 
#feas.dat<-read.csv("data/Suitability_v13_19.csv") #v13_19 updated with Ecologist review as of May 15, 2025
#feas.dat<-read.csv("data/Suitability_v13_22.csv") #v13_22 updated with Craig Delong review & inputed ratings Oct 1, 2025
feas.dat<-read.csv("data/Suitability_v13_24.csv") #v13_24 updated with Craig Delong & other review & non reviewed inputed ratings removed Dec 4, 2025

feas.dat<-subset(feas.dat, spp!='X')#remove any with no species defined 

feas.dat<-mutate(feas.dat, newsuit=if_else(is.na(newsuit), suitability, newsuit))#if not updated, use previous rating 


#take out the US and alberta stuff because it won't match plot data
feas.dat<-filter(feas.dat, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))

#load BEC plot data
load(file="data/tree_data_cleaned_updated.Rdata") 
#join BEC data with suitability data----
feas.dat.sub<-rename(feas.dat, ss_nospace_final=ss_nospace)
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
feas.dat.sub$Zone<-NULL #remove old column-incomplete 
feas.dat.sub<-mutate(feas.dat.sub, zone= case_when(grepl('ICH', ss_nospace)~"ICH",grepl('ESSF', ss_nospace)~"ESSF", grepl('MS', ss_nospace)~"MS", grepl('SWB', ss_nospace)~"SWB", 
                                                   grepl('SBPS', ss_nospace)~"SBPS", grepl('BAFA', ss_nospace)~"BAFA", grepl('CWH', ss_nospace)~"CWH", grepl('IDF', ss_nospace)~"IDF",
                                                   grepl('BG', ss_nospace)~"BG", grepl('ESSF', ss_nospace)~"ESSF", grepl('CDF', ss_nospace)~"CDF", grepl('SBS', ss_nospace)~"SBS", 
                                                   grepl('MH', ss_nospace)~"MH", grepl('CMA', ss_nospace)~"CMA", grepl('PP', ss_nospace)~"PP", grepl('BWBS', ss_nospace)~"BWBS", TRUE~ NA))

#calculate average abundances by feas scores
avgs<-group_by(feas.dat.sub, zone, bgc, ss_nospace, Species, spp, newsuit_ord)%>%
  summarise(mean_abund_ss=mean(TotalAB, na.rm = T), sd_abund_ss=sd(TotalAB, na.rm = T), nplots_ss=n())
avgs<-mutate(avgs, sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))

#spp plots ----
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

#QA/QC----
#filter plot data for quality 
#throw out anything with only 1 plot or high variability but low coverage (i.e. mean<sd & <5 plots)?
avgs<-mutate(avgs, sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))
avgs<-subset(avgs, nplots_ss>1) 
avgs$diff<-avgs$mean_abund_ss-avgs$sd_abund_ss
avgs<-mutate(avgs, remove=if_else(nplots_ss<5 & diff<0, 'Y', 'N')) 
avgs<-subset(avgs, remove=='N') 
avgs$remove<-NULL
avgs$diff<-NULL

#look at prop of site series out of total by spp 
nss<-select(feas.dat, spp)%>% group_by(spp)%>%summarise(n_site_series=n())
avgs<-group_by(avgs, spp)%>%mutate(n_ss=n())%>%left_join(., nss)%>%mutate(prop_ss= n_ss/n_site_series)
select(avgs, spp, prop_ss)%>%distinct(.)

max(avgs$nplots_ss)
mean(avgs$nplots_ss)
median(avgs$nplots_ss)

#pull out ss to review 
#check ratings with relative cutoffs
avgs<-mutate(avgs,review= case_when(newsuit_ord=="5" & mean_abund_ss>0~"Y",
                                    newsuit_ord=="4" & mean_abund_ss>1~"Y",
                                    newsuit_ord=="3" & mean_abund_ss> 10~"Y",
                                    newsuit_ord=="3" & mean_abund_ss< 1~"Y",
                                    newsuit_ord=="2" & mean_abund_ss> 25~"Y",
                                    newsuit_ord=="2" & mean_abund_ss< 10~"Y",
                                    newsuit_ord=="1" & mean_abund_ss<25~"Y",
                                    TRUE~"N"))
check<-subset(avgs, review=="Y")

feas.dat.sub<-left_join(feas.dat.sub, avgs)
feas.dat.validate<-subset(feas.dat.sub, !is.na(review))

#don't review parkland subzones, should all be E3/E4
feas.dat.validate<- mutate(feas.dat.validate, review= ifelse(grepl("p$", bgc), "N", review))


#create edatopes----
#create a unique variable for the 15 edatopic spaces for modeling 
feas.dat.validate<- mutate(feas.dat.validate, NutrientRegime_clean=ifelse(NutrientRegime_clean=="F","E", NutrientRegime_clean))# there is no F
feas.dat.validate<- mutate(feas.dat.validate, MoistureRegime_clean=ifelse(MoistureRegime_clean=="8","7", MoistureRegime_clean)) #lump 8s with 7s to match CCISS

feas.dat.validate$edatope<-paste(feas.dat.validate$NutrientRegime_clean, feas.dat.validate$MoistureRegime_clean, sep="")
feas.dat.validate<- mutate(feas.dat.validate, edatopex=case_when(edatope=="C3"|edatope=="C4"~"C34",
                                                                 edatope=="C1"|edatope=="C2"|edatope=="C0"~"C12",
                                                                 edatope=="C5"|edatope=="C6"|edatope=="C7"~"C56",
                                                                 edatope=="A3"|edatope=="A4"|edatope=="B3"|edatope=="B4"~"AB34",
                                                                 edatope=="A0"|edatope=="B0"|edatope=="A1"|edatope=="A2"|edatope=="B1"|edatope=="B2"~"AB12",
                                                                 edatope=="A7"|edatope=="B7"|edatope=="A5"|edatope=="A6"|edatope=="B5"|edatope=="B6"~"AB56",
                                                                 edatope=="D3"|edatope=="D4"|edatope=="E3"|edatope=="E4"~"DE34",
                                                                 edatope=="D0"|edatope=="E0"|edatope=="D1"|edatope=="D2"|edatope=="E1"|edatope=="E2"~"DE12",
                                                                 edatope=="D7"| edatope=="E7"|edatope=="D5"|edatope=="D6"|edatope=="E5"|edatope=="E6"~"DE56",
                                                                 TRUE~NA))
feas.dat.validate$edatope<-as.factor(feas.dat.validate$edatope)
sort(unique(feas.dat.validate$edatope))
sort(unique(feas.dat.validate$edatopex))

#remove unused cols 
names(feas.dat.validate)
feas.dat.validate<-select(feas.dat.validate, -X, -SuccessionalStatus, -PlotRepresenting, -SiteUnit, -elev_check, -BECSiteUnit)


#save output----
save(feas.dat.validate, file="data/feas_abund_data_validate.Rdata")


#set moist and nutrients as ordinal
#feas.dat.sub$MoistureRegime_clean<-ordered(feas.dat.sub$MoistureRegime_clean, levels = c(8,7,6,5,4,3,2,1,0))
#feas.dat.sub$NutrientRegime_clean<-ordered(feas.dat.sub$NutrientRegime_clean, levels = c("F", "E", "D", "C", "B", "A"))

#Look at deviations between expert ratings and plot data 
#-By edaphic space
#-By expert 
#-By species 
#-By expected vs actual deviation

ggplot(check, aes(x = newsuit_ord, y = mean_abund_ss, alpha=0.5))+
  #geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ spp, scales='free_y') + theme_bw() + theme(legend.position='none') 

ggplot(feas.dat.validate, aes(x = newsuit_ord, y = TotalAB, alpha=0.5))+
  #geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ edatopex) + theme_bw() + theme(legend.position='none') 


ggplot(subset(feas.dat.validate,review=="Y"), aes(y = spp))+ geom_bar()
ggplot(subset(feas.dat.validate,review=="Y"), aes(y = zone))+ geom_bar()
ggplot(subset(feas.dat.validate,review=="Y"), aes(y = edatope))+ geom_bar()
ggplot(subset(feas.dat.validate,review=="Y"), aes(y = mod))+ geom_bar()


  
#pull in climate data ----
##select climate variables from BGC model
#climrVars = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
#              "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
#              "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
#              "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", "PPT_JAS", "CMD.total")
#names(feas.dat) 

#from climr_getdata_plots.R
#load("data/clim_dat.plots.Rdata")
#clim_dat<-rename(clim_dat, PlotNumber=id)%>% select(c("PlotNumber", "PERIOD"), climrVars)

#feas.dat.clim<-left_join(feas.dat.sub, clim_dat)

#save(feas.dat.clim, file="data/feas_abund_clim_data.Rdata")

#calculate PC axes for climate params----
#source("scripts/Bayesian_analysis/climPCAs.R")



