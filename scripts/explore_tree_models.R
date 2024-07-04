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

#load QA/QCed tree data
rm(list = ls())
load(file="data/tree_data_cleaned.Rdata") #all
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
hist(Cw$TotalA_scaled)

#model 
library(brms)
library(tidyverse)

#model formula
#start simple with PPT and MAT only 
modform <- bf(TotalA_scaled~ MAT + PPT + (MAT + PPT||MoistureRegime_clean) +  (MAT + PPT||NutrientRegime_clean))
              

#set priors
priors <- c(prior(normal(0,1),class=b),
            prior(normal(0,1),class=b, coef="MAT"),
            prior(normal(0,1),class=b, coef="PPT"),
            prior(normal(0,1),class=Intercept),
            prior(cauchy(0,0.5), class = sd), 
            prior(cauchy(0,0.5), class = phi))

set_prior("normal(0,1)", class = "b")

#run in brms for each spp  
mod.Lw <- brm(modform, Lw ,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
                  control = list(adapt_delta=0.9, max_treedepth = 11), 
                  iter=6000, warmup = 1000, family = Beta(), inits = 0)
                  #file=".rmd") 
summary(mod.Cw)
MCMCvis::MCMCtrace(mod.Cw)


pp_check(mod.Cw) #looks good 


mod<-lmer(TotalA~ (1|bgc:SiteUnit) + (1|bgc) + (1|year) + (1|PlotNumber), Cw)
ranef(mod)
library(brms)

TotalA~ scale(MAT) + scale(PPT) +
  Aspect + Slope + (|GIS_BGC) + (MAT|GIS_BGC:SiteUnit),  Hw) #can't do nested RE

#categorize MAT?           
Hw<-mutate(Hw, MATcat=cut(MAT, breaks = 5,
                          labels = c("MAT1", "MAT2", "MAT3", "MAT4", "MAT5")))
Hw<-left_join(Hw, spp_tab)%>%subset(nobs>2)

Hwmod<-lmer(log(TotalA)~ scale(MAT) + scale(PPT) + (MATcat|bgc:Site) ,  Hw) #can't do nested Random slope 
summary(Hwmod)
mcmc_plot(mod.Cw)
Hwmod<-lmer(log(TotalA)~ scale(MAT) + scale(PPT) + (MATcat|GIS_BGC) + (MATcat|GIS_BGC:SiteSeries),  Hw) #can't do nested RE

#climate, SNR, aSMR, aspect, slope 

