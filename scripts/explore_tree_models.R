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

#load QA/QCed tree data----
rm(list = ls())
load(file="data/tree_data_cleaned.Rdata") #all
gc()

#filter to 16
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
spp_keep<-spp_keep$Species #16 
tree_dat<-subset(tree_dat, Species %in% spp_keep)
unique(tree_dat$Species)              

#load(file="data/tree_spp_data_cleaned.Rdata") #top 16 spp individually

#look at response variable
hist(tree_dat$TotalA)

#try beta distribution for total A response because it is proportional data (1-100)
#scale these to 0,1 so can use Beta dist 
tree_dat$TotalA_scaled<-(tree_dat$TotalA)/100

At$TotalA_scaled<-(At$TotalA)/100
Ba$TotalA_scaled<-(Ba$TotalA)/100
Bl$TotalA_scaled<-(Bl$TotalA)/100
Cw$TotalA_scaled<-(Cw$TotalA)/100
Ep$TotalA_scaled<-(Ep$TotalA)/100
Fd$TotalA_scaled<-(Fd$TotalA)/100
Hm$TotalA_scaled<-(Hm$TotalA)/100
Hw$TotalA_scaled<-(Hw$TotalA)/100
Lw$TotalA_scaled<-(Lw$TotalA)/100
Pl$TotalA_scaled<-(Pl$TotalA)/100
Py$TotalA_scaled<-(Py$TotalA)/100
Sb$TotalA_scaled<-(Sb$TotalA)/100
Se$TotalA_scaled<-(Se$TotalA)/100
Ss$TotalA_scaled<-(Ss$TotalA)/100
Sw$TotalA_scaled<-(Sw$TotalA)/100
Yc$TotalA_scaled<-(Yc$TotalA)/100

#check
hist(tree_dat$TotalA_scaled)

#can also try a normal distribution w/ logged response b/c easier to interpret results- but fails on homoscedasticity of residuals 
hist(log(tree_dat$TotalA))
tree_dat$TotalA_log<-log(tree_dat$TotalA)

#random forest models----
#get all climate variables 
vars<-climr::variables #look up table for vars 
var_names<-vars$Code
#add response variable and a few other preds 
var_names<-c("TotalA", "Species", "PlotNumber",  "NutrientRegime_clean", "MoistureRegime_clean", "SlopeGradient", "Aspect", "year", "bgc" , var_names)

#create dataset with all preds for rf model 
tree_dat_sub<-select(tree_dat, var_names)
tree_dat_sub<-na.omit(tree_dat_sub)  #remove NAs

#run model
RFmodel_allspp <- ranger::ranger(
  TotalA~ .,
  data = tree_dat_sub,
  splitrule = "maxstat",  
  importance = 'permutation',
  scale.permutation.importance = TRUE, 
  mtry=253)
print(RFmodel_allspp) #R2=0.22

#look at importance of features
importancedf<-as.data.frame(importance(RFmodel_allspp))

#plot climate variables----
#plot by edatopic grid space 
#look at factors marked important in RF
# RH (whole year)
#summer PPT
#Eref autumn 
#TD
#PAS (whole year)- less impt

ggplot(tree_dat, aes(y=TotalA, x=Eref_at, fill=MoistureRegime_clean))+ 
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  #geom_smooth()  +
  ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=PPT_sm))+ 
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth()  +
  ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=RH))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)


ggplot(tree_dat, aes(y=TotalA, x=TD))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=SlopeGradient, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=PAS, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=year, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=as.factor(year)))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=bgc))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

#check for collinearity among climate and enviro preds
check_climate<-select(tree_dat, Eref_at, PAS, RH, PPT_sm, TD, SlopeGradient)
pairs(check_climate)
#RH and TD look highly correlated

ggplot(tree_dat, aes(y=RH, x=TD))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  
cor.test(tree_dat$RH, tree_dat$TD)#86% overall - too high- just use RH since higher importance

#bayesian models---- 
library(brms)
library(tidyverse)

#mean center continuous vars
#tree_dat$DD5_scaled<-scale(tree_dat$DD5)
#hist(tree_dat$DD5_scaled)
#tree_dat$CMI_scaled<-scale(tree_dat$CMI)
#hist(tree_dat$CMI_scaled)
#tree_dat$TD_scaled<-scale(tree_dat$TD)
#hist(tree_dat$TD_scaled)
#tree_dat$PPT_10_scaled<-scale(tree_dat$PPT_10)
#hist(tree_dat$PPT_10_scaled)
tree_dat$SlopeGradient_scaled<-scale(tree_dat$SlopeGradient)
hist(tree_dat$SlopeGradient_scaled)
#tree_dat$Aspect_scaled<-scale(tree_dat$Aspect)
#hist(tree_dat$Aspect_scaled)
tree_dat$RH_scaled<-scale(tree_dat$RH)
tree_dat$Eref_at_scaled<-scale(tree_dat$Eref_at)
tree_dat$PPT_sm_scaled<-scale(tree_dat$PPT_sm)
tree_dat$PAS_scaled<-scale(tree_dat$PAS)
tree_dat$year_scaled<- scale(tree_dat$year) #for changes over time in plot cover 

#make a factor for year as well 
tree_dat$year_factor<- as.factor(tree_dat$year) #for interannual differences in sampling intensity etc. - phi model 

#set priors----
#use set prior option to set to all coefs 
beta_priors <- c(set_prior("normal(0,1)", class = "Intercept"),  
            set_prior("normal(0, 0.5)", class = "b"), 
            set_prior("cauchy(0,0.5)", class = "sd"))

priors <- c(set_prior("normal(2,1)", class = "Intercept"), #on logged intercept 
            set_prior("normal(0, 0.5)", class = "b"), 
            set_prior("cauchy(0,0.5)", class = "sd"), 
            set_prior("cauchy(0, 0.5)", class = "sigma"))


#model formulas----
#v0- species separately 
modform0<-bf(TotalA_log~ TD_scaled + PPT_10_scaled + SlopeGradient_scaled + Aspect_scaled + (TD_scaled + PPT_10_scaled + SlopeGradient_scaled + Aspect_scaled ||MoistureRegime_clean) +  
     (TD_scaled + PPT_10_scaled + SlopeGradient_scaled + Aspect_scaled ||NutrientRegime_clean))

# v1- using vars suggested by Kiri/Will DD5, CMI
#beta dist
modform <- bf(TotalA_scaled~ DD5_scaled + CMI_scaled + (DD5_scaled + CMI_scaled||Species:MoistureRegime_clean) +  
                                        (DD5_scaled + CMI_scaled||Species:NutrientRegime_clean))

# v2-try predictors found strongest Wang et al. 2012 Random Forest= TD, PPT_10
#https://www2.gov.bc.ca/assets/gov/environment/natural-resource-stewardship/nrs-climate-change/applied-science/wangfinalreport.pdf
#beta dist
modform2 <- bf(TotalA_scaled~ TD_scaled + PPT_10_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))
# v3- normal dist, same preds
modform3 <- bf(TotalA_log~ TD_scaled + PPT_10_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))
# v4- normal dist, add slope and aspect 
modform4 <- bf(TotalA_log~ TD_scaled + PPT_10_scaled + SlopeGradient_scaled + Aspect_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))
# v5- normal dist- RH, Eref_at, PPT_sm, PAS - found by importance ranking from RF. plus year
#should bgc go in this model?? if so, in nested group term or main model or both?? or nowhere bc climate already informs??
modform5 <- bf(TotalA_log~ Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled + SlopeGradient_scaled + year_scaled + (Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled||Species:MoistureRegime_clean) +  
                 (Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled||Species:NutrientRegime_clean))
# v6- same as v5 but beta dist plus year_factor + species in phi model and remove year from main model
modform6 <- bf(TotalA_scaled~ Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled + (Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled||Species:MoistureRegime_clean) +  
                 (Eref_at_scaled + PPT_sm_scaled + RH_scaled + PAS_scaled||Species:NutrientRegime_clean),
                 phi~ Species + year_factor)

#run in brms----
#all species together 
#update model number and file name when running diff versions 
mod.all6 <- brm(modform6, tree_dat ,cores=6, chains=3, threads = threading(12), backend = "cmdstanr", prior = beta_priors, 
                  #control = list(adapt_delta=0.99, max_treedepth = 11), 
                  iter=6000, warmup = 1000, init = 0, family = Beta(),
                file= "outputs/brms/mod_allspp6.Rmd") 
summary(mod.all5)


#try species separately??
Hw<-subset(tree_dat, Species=='TSUGHET') #4716 obs 
mod.HW<-brm(modform0, Hw ,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
            #control = list(adapt_delta=0.99, max_treedepth = 11), 
            iter=6000, warmup = 1000, init = 0, 
            file= "outputs/brms/mod_Hw.Rmd")
