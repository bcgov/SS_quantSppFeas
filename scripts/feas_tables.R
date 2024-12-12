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
#pulls in cleaned feasibility tables w/ climate data created from 0_DataPrep_OrdinalForests.Rmd in Tree Feasibility Prediction repr @whmacken 
#filters BGCs further
#runs PCA on climate vars and merges with feas values 
#merges feas values with BEC plot data (tree abundance)
#calculates and plots mean abundances by feas rating for top spp 

#load feasibility data with site & climate info---- 
feas.dat <- readRDS("C:/Users/ccollins/OneDrive - Government of BC/WNA_BGC/OrdinalForest_data.rds")

#filter out small BGCs not included in training/projections (from BGC projections repo- create training set.R)
BGC_exclude<-c( "BWBScm"  ,"CMAun_OR",  "CMAwh",   "CWHxs",  "ESSFdcp", "ESSFdh1" ,"ESSFdh2"  ,"ESSFdh_WA","ESSFwh_MT","ESSFxcp", "ESSFxcw",  
                "ESSFxh_WA", "ESSFxvw" ,  "ICHmc1a",  "ICHxwa"  ,  "IDFww"  ,   "IDFww1"   , "IDFww2"  ,  "IDFxx1"  ,  "MHun"  ,    "MHunp"  ,   "MHwh" ,    
                "MHwhp"  ,   "MSdc2"  ,   "MSun")

feas.dat<-subset(feas.dat, !(BGC %in% BGC_exclude))
feas.dat<-na.omit(feas.dat)

unique(feas.dat$BGC)#309-> 366 in projections, missing 57 subzones 
unique(feas.dat$spp)
unique(tree_dat$Species)

#select climate variables from BGC model
climrVars = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
              "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
              "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
              "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", "PPT_JAS", "CMD.total")
names(feas.dat)

feas.dat<-mutate(feas.dat, CMD.total=CMD.def +CMD)
varsl<-c(c("zone", "BGC", "SS_NoSpace","spp", "newfeas",  "WNA_DEM_4326_clipped" ,"xcoord", "ycoord") , climrVars)
feas.dat.sub<-select(feas.dat, varsl)

#calculate PC axes for climate params----
library(vegan)
library(ape)

climdat<- select(feas.dat.sub, climrVars) %>%
  #mutate(PASlog=log(PAS), PAS_splog=log(PAS_sp),PPT_MJlog=log(PPT_MJ), PPT_JASlog=log(PPT_JAS))%>% #log precip vars
  #select(-PAS, -PAS_sp, -PPT_MJ, -PPT_JAS) %>%
  scale(.) # standardize the data

# Perform principal component analysis (PCA) on the standardized data
climpca<-  prcomp(climdat)

# Print a summary of the PCA results
summary(climpca)
biplot(climpca)

# Convert the PCA results to a tibble and add a column for observation names
climpca_df <- as_tibble(climpca$x) %>% 
  mutate(observation = rownames(climpca$x))%>%select(PC1, PC2, PC3, observation)

#in vegan
climpca2 <- vegan::rda(climdat, scale=F) #already scaled above 
scores <- scores(climpca)
biplot(climpca) 
eigenvals(climpca2)
#loading scores of climate params onto PCA axes 
climloads<-envfit(climpca, climdat, permutations = 999, strata = NULL)
climloads
scrs <- as.data.frame(scores(climloads, display = "vectors"))
scrs<-add_rownames(scrs, "clim")

ggplot(scrs) +
  geom_point(mapping = aes(x = PC1, y = PC2)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = scrs,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = scrs, aes(x = PC1, y = PC2, label = clim),
            size = 3)

#merge PCs back with feas data---- 
feas.dat.sub<-add_rownames(feas.dat.sub, "observation")
feas.dat.sub<-left_join(climpca_df, feas.dat.sub)

save(feas.dat.sub, file="data/feasibility_data.Rdata")


#merge feas & tree data---- 
feas.dat.sub<-mutate(feas.dat.sub, Species= case_when(spp=="Ba"~"ABIEAMA",
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
                                                      spp=="Plc"~"PINUCON", #coastal
                                                      spp=="Pli"~"PINUCON", #interior
                                                      spp=="Pyc"~"PINUPON",
                                                      spp=="Pyi"~"PINUPON",
                                                      spp=="Pw"~"PINUMON",
                                                      spp=="Acb"~"POPUBAL",
                                                      spp=="At"~"POPUTRE",
                                                      spp=="Act"~"POPUTRI",
                                                      spp=="Fd"~"PSEUMEN", 
                                                      spp=="Fdi"~"PSEUMEN", #interior
                                                      spp=="Fdc"~"PSEUMEN", #coastal
                                                      spp=="Tw"~ "TAXUBRE",
                                                      spp=="Cw"~ "THUJPLI",
                                                      spp=="Hw"~ "TSUGHET",
                                                      spp=="Hm"~"TSUGMER")) 
tree_dat_sub<-select(tree_dat, PlotNumber, Species, TotalA, SiteUnit, bgc,NutrientRegime_clean,MoistureRegime_clean, Latitude, Longitude)%>%rename(BGC=bgc)

tree_dat_sub<-mutate(tree_dat_sub,  SS_NoSpace= gsub(" ", "", SiteUnit)) #create matching column to feas table

names(feas.dat.sub)
feas.dat.sub<-left_join(feas.dat.sub, tree_dat_sub)  
#check merge worked ok 
plot(feas.dat.sub$ycoord~feas.dat.sub$Latitude)??
  
avgs<-group_by(feas.dat.sub, zone, BGC, SS_NoSpace, Species, spp, newfeas_ord)%>%summarise(mean_abund_ss=mean(TotalA, na.rm = T))%>%ungroup(.)%>%
group_by(BGC)%>% mutate(mean_abund_subzone=mean(mean_abund_ss, na.rm = T))%>% ungroup(.)%>%
group_by(zone)%>% mutate(mean_abund_zone=mean(mean_abund_subzone, na.rm = T))#4 is actually a 5 here (i.e. zero)

avgs<-na.omit(avgs)



#plot by spp with data 
check<-group_by(feas.dat.sub, spp)%>%summarise(n_ratings= n_distinct(newfeas), abund=mean(TotalA, na.rm=T))
check<-subset(check, n_ratings>1 & !is.na(abund))
check$spp

#At- 1 vs 2 not clear in IDF
ggplot(subset(avgs, Species=="POPUTRE"), aes(x = newfeas_ord, y = log(mean_abund_ss+1), colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Ba- fine, rare
ggplot(subset(avgs, Species=="ABIEAMA"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Bl-looks good
ggplot(subset(avgs, Species=="ABIELAS"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Cw- looks good 
ggplot(subset(avgs, Species=="TSUGHET"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Ep- looks good 
ggplot(subset(avgs, Species=="BETUPAP"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Fd - 4s too high in most zones - driven by Fdi & Fdc entries 
ggplot(subset(avgs, Species=="PSEUMEN"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
ggplot(subset(avgs, Species=="PSEUMEN"& spp=="Fd"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
ggplot(subset(avgs, spp=="Fdi"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
ggplot(subset(avgs, spp=="Fdc"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Hm- fine, rare 
ggplot(subset(avgs, Species=="TSUGMER"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Hw- looks good 
ggplot(subset(avgs, Species=="TSUGHET"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Lw-looks good 
ggplot(subset(avgs, Species=="LARIOCC"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Pl- 4s too high in most zones- driven by Plc and Pli entries... 
ggplot(subset(avgs, Species=="PINUCON"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Plc, Plis have high abundances for 4s 
ggplot(subset(avgs, spp=="Plc"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() #+ theme(legend.position='none')
ggplot(subset(avgs, spp=="Pli"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() #+ theme(legend.position='none')
#remove them 
ggplot(subset(avgs, Species=="PINUCON"& spp=="Pl"), aes(x = newfeas_ord, y = mean_abund_ss, alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Ss- looks good 
ggplot(subset(avgs, Species=="PICEMAR"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Ss- high values for 4 in CWH
ggplot(subset(avgs, Species=="PICESIT"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')
#Yc - 1 & 2 almost the same in CWH
ggplot(subset(avgs, Species=="CALLNOO"), aes(x = newfeas_ord, y = mean_abund_ss, colour="blue", alpha=0.5))+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white", position=position_dodge(width=0.9), alpha=0.5) +
  facet_wrap( ~ zone) + theme_bw() + theme(legend.position='none')

##OLD###
#filter feas table
#take out the US and alberta stuff because it won't match plot data
feas_tab<-filter(feas_tab, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abE|abC|SBAP|SASbo|PPxh|MSd|MSx', bgc))
#subset to only top 16 spp
spp_tab0<-tree_dat%>%  group_by(Species)%>%  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
feas_tab<-subset(feas_tab, Species %in% spp_keep$Species)
rm(spp_tab0)


spp_keep$Species


