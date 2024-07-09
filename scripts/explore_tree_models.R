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
#filter to 16
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
spp_keep<-spp_keep$Species #16 
tree_dat<-subset(tree_dat, Species %in% spp_keep)
unique(tree_dat$Species)              

load(file="data/tree_spp_data_cleaned.Rdata") #top 16 spp individually

#look at response variable
hist(tree_dat$TotalA)

hist(Cw$TotalA)

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

#can also try a normal distribution w/ logged response b/c easier to interpret results
hist(log(tree_dat$TotalA))
tree_dat$TotalA_log<-log(tree_dat$TotalA)

#plot climate variables----
#plot by edatopic grid space 
ggplot(tree_dat, aes(y=TotalA, x=MAT, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth()  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=PPT, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)
#remove outliers for PINUPON and PICEMAR 
tree_dat<-subset(tree_dat, ID!=-1315162869 & ID!=1793185075)

ggplot(tree_dat, aes(y=TotalA, x=AHM, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=CMD, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=CMI, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)

ggplot(tree_dat, aes(y=TotalA, x=DD5, fill=MoistureRegime_clean))+
  geom_point() + facet_wrap(~Species, scales = 'free') + 
  geom_smooth(method='lm')  + ylim(0,100)


#model---- 
library(brms)
library(tidyverse)

#mean center continuous vars
tree_dat$DD5_scaled<-scale(tree_dat$DD5)
hist(tree_dat$DD5_scaled)
tree_dat$CMI_scaled<-scale(tree_dat$CMI)
hist(tree_dat$CMI_scaled)
tree_dat$TD_scaled<-scale(tree_dat$TD)
hist(tree_dat$TD_scaled)
tree_dat$PPT_10_scaled<-scale(tree_dat$PPT_10)
hist(tree_dat$PPT_10_scaled)
tree_dat$SlopeGradient_scaled<-scale(tree_dat$SlopeGradient)
hist(tree_dat$SlopeGradient_scaled)
tree_dat$Aspect_scaled<-scale(tree_dat$Aspect)
hist(tree_dat$Aspect_scaled)

#model formulas----
# v1- using vars suggested by Kiri/Will DD5, CMI
#beta dist
#modform <- bf(TotalA_scaled~ DD5_scaled + CMI_scaled + (DD5_scaled + CMI_scaled||Species:MoistureRegime_clean) +  
#                                        (DD5_scaled + CMI_scaled||Species:NutrientRegime_clean))

# v2-try predictors found strongest Wang et al. 2012 Random Forest= TD, PPT_10
#https://www2.gov.bc.ca/assets/gov/environment/natural-resource-stewardship/nrs-climate-change/applied-science/wangfinalreport.pdf
#beta dist
#modform <- bf(TotalA_scaled~ TD_scaled + PPT_10_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
#                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))
# v3- normal dist
#modform <- bf(TotalA_log~ TD_scaled + PPT_10_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
#                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))
# v4- try with slope and aspect 
modform <- bf(TotalA_log~ TD_scaled + PPT_10_scaled + SlopeGradient_scaled + Aspect_scaled + (TD_scaled + PPT_10_scaled||Species:MoistureRegime_clean) +  
                (TD_scaled + PPT_10_scaled||Species:NutrientRegime_clean))

#set priors----
#priors <- c(prior(normal(0,1),class=b),
            #prior(normal(0,1),class=b, coef="DD5_scaled"),
            #prior(normal(0,1),class=b, coef="CMI_scaled"),
            #prior(normal(0,1),class=b, coef="TD_scaled"),
            #prior(normal(0,1),class=b, coef="PPT_10_scaled"),
#            prior(normal(0,1),class=Intercept),
#            prior(cauchy(0,0.5), class = sd))#, 
            #prior(cauchy(0,0.5), class = phi))

#use set prior option to set to all coefs 
priors <- c(set_prior("normal(0,1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"),
            set_prior("cauchy(0,0.5)", class = "sd"), 
            set_prior("normal(0,1)", class = "sigma"))

#run in brms----
#all species together 
#update model number and file name when running diff versions 
mod.all4 <- brm(modform, tree_dat ,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
                  #control = list(adapt_delta=0.99, max_treedepth = 11), 
                  iter=6000, warmup = 1000, init = 0, 
                  file= "mod_allspp4.Rmd") 
summary(mod.all2)
MCMCvis::MCMCtrace(mod.all)


pp_check(mod.all2) #looks okay 



#climate, SNR, aSMR, aspect, slope 
