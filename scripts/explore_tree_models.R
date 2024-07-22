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

#libraries 
library(tidyverse)

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

#check
hist(tree_dat$TotalA_scaled)

#can also try a normal distribution w/ logged response b/c easier to interpret results- but fails on homoscedasticity of residuals 
hist(log(tree_dat$TotalA))
tree_dat$TotalA_log<-log(tree_dat$TotalA)

#plot climate variables----

#get all climate variables 
vars<-climr::variables #look up table for vars 
var_names<-vars$Code
#add response variable and a few other preds to all climate data 
#var_names<-c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "SlopeGradient", "Aspect", "year", "bgc" , "Elevation" ,var_names)

#subset predictor vars
var_names<-c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "year", "bgc" , "year", "bgc","Tmax_sm", "TD", "PPT_sm", "DD5_sp")

#create dataset with all preds for rf model 
tree_dat_sub<-select(tree_dat, var_names)
tree_dat_sub<-na.omit(tree_dat_sub)  #remove NAs

#check for correlation among climate vars
pairs(tree_dat_sub[,c(7:10)])

plot(tree_dat$DD5, tree_dat$Tmax_sm)
cor.test(tree_dat$DD5, tree_dat$Tmax_sm)#correlated ~0.7 remove
cor.test(tree_dat$DD5_sp, tree_dat$Tmax_sm)#correlated ~0.65 better
plot(tree_dat$CMI, tree_dat$TD)
cor.test(tree_dat$NFFD, tree_dat$TD) #r~ 0.45 ok 

#random forest models----
#run model
RFmodel_allspp <- ranger::ranger(
  TotalA~ .,
  data = tree_dat_sub,
  splitrule = "maxstat",  #recommended for accurate importance feature ranking
  importance = 'permutation',
  mtry=11, 
  scale.permutation.importance = TRUE)

save(RFmodel_allspp, file="outputs/ranger/RFmodel5.Rdata")
print(RFmodel_allspp) #R2 max=0.3

#look at importance of features
importancedf<-as.data.frame(ranger::importance(RFmodel_allspp))

# Split the dataset into training and testing sets
set.seed(123) 
train_indices <- sample(1:nrow(tree_dat_sub), 0.7 * nrow(tree_dat_sub)) # 70% of the data for training
train_data <- tree_dat_sub[train_indices, ]
test_data <- tree_dat_sub[-train_indices, ]

#try running again on training data to get more accurate R2 value
RFmodel_allspp <- ranger::ranger(
  TotalA~ .,
  data = train_data,
  splitrule = "maxstat",  #recommended for accurate importance feature ranking
  importance = 'permutation',
  mtry=9, 
  scale.permutation.importance = TRUE)

save(RFmodel_allspp, file="outputs/ranger/RFmodel7.Rdata")
print(RFmodel_allspp) #R2 max=0.3

# Predict on the test set
preds <- predict(RFmodel_allspp, data = test_data)

#regress
plot(test_data$TotalA, preds$predictions)
cor.test(test_data$TotalA, preds$predictions)#r=0.52, r2=0.27 (so about the same as OOB)


#ordinal forest models----
# Convert totalA to an ordered factor
#breaks from plot data (explore_veg_data.R)
#newfeas  meds  mean
#       1    17 21.8 
#       2    10 16.2 
#       3     5  9.64
#       4     2  2   
#       5     2  4.67

tree_dat_sub<-mutate(tree_dat_sub, cover_rank=case_when(TotalA>20~1, 
                                                         TotalA<7.5~3, 
                                                         TRUE~2))
hist(tree_dat_sub$cover_rank)

group_by(tree_dat_sub, cover_rank)%>%summarise(counts=n())#about even

#set as ordinal factor
tree_dat_sub$cover_rank<-ordered(tree_dat_sub$cover_rank, levels = c(3, 2, 1))
str(tree_dat_sub)#good 

#remove continuous response
tree_dat_sub$TotalA<-NULL

# Split the dataset into training and testing sets
set.seed(123) 
train_indices <- sample(1:nrow(tree_dat_sub), 0.7 * nrow(tree_dat_sub)) # 70% of the data for training
train_data <- tree_dat_sub[train_indices, ]
test_data <- tree_dat_sub[-train_indices, ]

#Fit the ordinal random forest model
library(ordinalForest)
RFord <- ordfor(depvar = "cover_rank", data = train_data, mtry = 9, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodel2.Rdata")

# Predict on the test set
preds <- predict(RFord, newdata=test_data)

#confusion matrix 
conf_mat<-table('true'=test_data$cover_rank, 'predicted'=preds$ypred) #predicts data wrong about 50% of the time 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.5



#bayesian models---- 
library(brms)
library(tidyverse)

#mean center continuous vars
tree_dat$DD5_sp_scaled<-scale(tree_dat$DD5_sp)
hist(tree_dat$DD5_sp_scaled)
#tree_dat$CMI_scaled<-scale(tree_dat$CMI)
#hist(tree_dat$CMI_scaled)
tree_dat$TD_scaled<-scale(tree_dat$TD)
hist(tree_dat$TD_scaled)
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
tree_dat$Tmax_sm_scaled<-scale(tree_dat$Tmax_sm)

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
summary(mod_allspp6.Rmd)


#try species separately??
Hw<-subset(tree_dat, Species=='TSUGHET') #4716 obs 
mod.HW<-brm(modform0, Hw ,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
            #control = list(adapt_delta=0.99, max_treedepth = 11), 
            iter=6000, warmup = 1000, init = 0, 
            file= "outputs/brms/mod_Hw.Rmd")
pp_check(mod_allspp6.Rmd)


#try with ordinal response---- 
#https://bookdown.org/content/3686/ordinal-predicted-variable.html
#https://kevinstadler.github.io/notes/bayesian-ordinal-regression-with-random-effects-using-brms/

#transform response 
tree_dat<-mutate(tree_dat, cover_rank=case_when(TotalA>20~1, 
                                                        TotalA<7.5~3, 
                                                        TRUE~2))
#set as ordinal factor
tree_dat$cover_rank<-ordered(tree_dat$cover_rank, levels = c(3, 2, 1))

#set priors 
priors <- c(set_prior("normal(0,4)", class = "Intercept"), 
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 0.5)", class = "b"))
                     
modform_ord<-bf(cover_rank~ PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled+ (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled||Species:MoistureRegime_clean) +  
                  (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled||Species:NutrientRegime_clean))

mod_ord<-brm(modform_ord, tree_dat,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
    control = list(adapt_delta=0.99, max_treedepth = 11), 
    iter=6000, warmup = 1000, init = 0, family=  cumulative("probit"),
    file= "outputs/brms/mod_ord")
