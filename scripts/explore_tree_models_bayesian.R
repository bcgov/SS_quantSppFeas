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

#load libs and data---- 
library(brms) 
library(ordinal)
library(lme4)
library(glmmTMB)
library(tidyverse)
library(remotes)
library(NBZIMM) #zero infl gaussian
library(bhsdtr2) #hierarchical ordinal bayesian  

#load feasibility plus abundance dataset with climate data 
load(file="data/feasibility_abundance_data.Rdata")
#currently using CMI and Tave_sm because they load most heavily on PC1 and PC2 but can also just use PC axes 

#set moist and nutrients as ordinal
#feas.dat.sub$MoistureRegime_clean<-newfeas_ord<-ordered(feas.dat.sub$MoistureRegime_clean, levels = c(8,7,6,5,4,3,2,1,0))
#feas.dat.sub$NutrientRegime_clean<-newfeas_ord<-ordered(feas.dat.sub$NutrientRegime_clean, levels = c("F", "E", "D", "C", "B", "A"))

#create a unique variable for the 15 edatopic spaces  
feas.dat.sub$edatopex<-paste(feas.dat.sub$NutrientRegime_clean, feas.dat.sub$MoistureRegime_clean, sep="")
feas.dat.sub<- mutate(feas.dat.sub, edatope=case_when(edatopex=="C3"|edatopex=="C4"~"C34",
                                                      edatopex=="C1"|edatopex=="C2"~"C12",
                                                      edatopex=="C5"|edatopex=="C6"~"C56",
                                                      edatopex=="C0"~"C0",
                                                      edatopex=="C7"~"C7",
                                                      edatopex=="A3"|edatopex=="A4"|edatopex=="B3"|edatopex=="B4"~"AB34",
                                                      edatopex=="A1"|edatopex=="A2"|edatopex=="B1"|edatopex=="B2"~"AB12",
                                                      edatopex=="A5"|edatopex=="A6"|edatopex=="B5"|edatopex=="B6"~"AB56",
                                                      edatopex=="A0"|edatopex=="B0"~"AB0",
                                                      edatopex=="A7"|edatopex=="B7"~"AB7",
                                                      edatopex=="A0"|edatopex=="B0"~"AB0",
                                                      edatopex=="D3"|edatopex=="D4"|edatopex=="E3"|edatopex=="E4"~"DE34",
                                                      edatopex=="D1"|edatopex=="D2"|edatopex=="E1"|edatopex=="E2"~"DE12",
                                                      edatopex=="D5"|edatopex=="D6"|edatopex=="E5"|edatopex=="E6"~"DE56",
                                                      edatopex=="D0"|edatopex=="E0"~"DE0",
                                                      edatopex=="D7"| edatopex=="D8"|edatopex=="E7"~"DE78",  #8 not in BYBEC
                                                      edatopex=="F5"|edatopex=="F6"~"F56", #F not in ByBEC
                                                      edatopex=="F3"~"F3", #F not in BYBEC
                                                      TRUE~NA))
feas.dat.sub$edatope<-as.factor(feas.dat.sub$edatope)


#build prior model----
#remove plot data
prior_df<-select(feas.dat.sub, -TotalA, -PlotNumber, -Latitude, -Longitude)%>%distinct(.)
prior_df$Tave_sms<-c(scale(prior_df$Tave_sm)) #use c() so doesn't change class type 
prior_df$CMIs<-c(scale(prior_df$CMI))#use c() so doesn't change class type
hist(prior_df$Tave_sms)
hist(prior_df$CMIs)

#simulate cover data from feas ratings 
knitr::kable(group_by(prior_df, newfeas_ord)%>%summarise(counts=n())) 

# Define the ranges for each class
#use proposed cutoffs from Mariotte et al. 2013a
ranges <- list(
  feas1 = c(12.01, 99),  #12% is dominance - high
  feas2 = c(2, 12 ),  #subordinates - med
  feas3 = c(0.01, 1.99 ) #remaining non zero -low   
  #true 4s -i.e. trace could be 0-0.5 or something??
)

# Simulate random data for each class within the given range
feas1_data <- runif(1864, min = ranges$feas1[1], max = ranges$feas1[2]) #points from kable cts
feas2_data <- runif(2540, min = ranges$feas2[1], max = ranges$feas2[2])
feas3_data <- runif(2135, min = ranges$feas3[1], max = ranges$feas3[2])

# Combine into a data frame
sim_data <- data.frame(
  value = c(feas1_data, feas2_data, feas3_data),
  newfeas = rep(c("1", "2", "3"), 
              times = c(length(feas1_data),length(feas2_data), length(feas3_data))))

prior_df$sim_abund <- NA

# Assign simulated data values based on matching class labels
for (feas in unique(prior_df$newfeas)) {
  # Find the corresponding simulated data for this class
  sim_values <- sim_data$value[sim_data$newfeas == feas]
    # Match the class rows in `existing_data` and assign the simulated values
  prior_df$sim_abund[prior_df$newfeas == feas] <- sim_values
} #ignore error 

#put in zeroes for 4s
prior_df<-mutate(prior_df, sim_abund= if_else(is.na(sim_abund), 0, sim_abund))

#run for one spp (Fd)
Fdfeas<-subset(prior_df, spp=="Fd")
str(Fdfeas)

#FIT MODEL 
#continuous response - hurdle lognormal 
#brms 
Fd_priormod<-brm(bf(sim_abund ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope),hu ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope)),
  data = Fdfeas,
  family = hurdle_lognormal(),
  chains = 3, iter = 2000, warmup = 1000)
save(Fd_priormod, file= "outputs/brms/Fd_priormod.Rdata")


#look at model output
summary(Fd_priormod)
ranef(Fd_priormod)
pp_check(Fd_priormod)
pred <- posterior_predict(Fd_priormod)
bayesplot::ppc_dens_overlay(y = log1p(Fdfeas$sim_abund), 
                            yrep = log1p(pred[1:10,]))

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Fdfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))

conditional_effects(Fd_priormod, effects = "CMIs", conditions = conditions,
                    re_formula = NULL) 

conditional_effects(Fd_priormod, effects = "Tave_sms", conditions = conditions,
                       re_formula = NULL) 

#posterior model---- 
#subset to one spp 
Fdabund<-subset(feas.dat.sub, spp=="Fd")
names(Fdabund)
Fdabund$Tave_sms<-c(scale(Fdabund$Tave_sm)) #use c() so doesn't change class type 
Fdabund$CMIs<-c(scale(Fdabund$CMI))#use c() so doesn't change class type

#FIT MODEL first with default (flat) priors
Fd_postmod_flat<-brm(bf(TotalA ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope),hu ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope)),
                 data = Fdabund, 
                 family = hurdle_lognormal(),
                 chains = 3, iter = 2000, warmup = 1000)
save(Fd_postmod_flat, file= "outputs/brms/Fd_postmod_flat.Rdata")

summary(Fd_postmod_flat)
pp_check(Fd_postmod_flat)
pred <- posterior_predict(Fd_postmod_flat)
bayesplot::ppc_dens_overlay(y = log1p(Fdabund$TotalA), 
                            yrep = log1p(pred[1:10,]))


#set informed priors from prior model
get_prior(Fd_postmod_flat)
summary(Fd_priormod)

priors<- c(set_prior("normal(1.34, 0.19)", class = "Intercept"), 
          set_prior("normal(0.76, 0.55)", class = "Intercept", dpar='hu'),
          set_prior("normal(0.52, 0.10)", class = "b", coef = "Tave_sms"),
          set_prior("normal(-0.66, 0.24)", class = "b", coef = "CMIs"),
          set_prior("normal(-1.81, 0.15)", class = "b", coef = "Tave_sms", dpar='hu'),
          set_prior("normal(0.84, 0.18)", class = "b", coef = "CMIs", dpar='hu'),
          set_prior("normal(0.47, 0.18)", class = "sd", group = "edatope", coef = "Intercept"), 
          set_prior("normal(0.14, 0.11)", class = "sd", group = "edatope", coef = "Tave_sms"),   
          set_prior("normal(0.48, 0.23)", class = "sd", group = "edatope", coef = "CMIs"),        
          set_prior("normal(2.10, 0.58)", class = "sd", group = "edatope", dpar='hu', coef = "Intercept"),  
          set_prior("normal(0.23, 0.17)", class = "sd", group = "edatope", dpar='hu', coef = "Tave_sms"),   
          set_prior("normal(0.30, 0.22)", class = "sd", group = "edatope", dpar='hu', coef = "CMIs"))      

#FIT MODEL again with expert informed priors
Fd_postmod_inf<-brm(bf(TotalA ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope),hu ~ Tave_sms + CMIs + (Tave_sms + CMIs||edatope)),
                     data = Fdabund, prior=priors,
                     family = hurdle_lognormal(),
                     chains = 3, iter = 2000, warmup = 1000)
save(Fd_postmod_inf, file= "outputs/brms/Fd_postmod_inf.Rdata")
get_prior(Fd_postmod_inf)

#look at model output
summary(Fd_postmod_inf)
ranef(Fd_postmod_inf)
pp_check(Fd_postmod_inf)
pred <- posterior_predict(Fd_postmod_inf)
bayesplot::ppc_dens_overlay(y = log1p(Fdabund$TotalA), 
                            yrep = log1p(pred[1:10,]))
#not much difference in pp_check compared to post-flat prior model :/ not sure if priors are actually getting set correctly or defaults still being used?.. 12/19/24
#slightly overpredicting at lower abundances and underpredicting at higher abundances... how to update priors or sim data to address this? 

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Fdabund$edatope)) |> 
  mutate(cond__ = paste0(edatope))

conditional_effects(Fd_postmod_inf, effects = "CMIs", conditions = conditions,
                    re_formula = NULL, prob = 0.9)
#compare to default priors -different!
conditional_effects(Fd_postmod_flat, effects = "CMIs", conditions = conditions,
                    re_formula = NULL, prob=0.9) 

conditional_effects(Fd_postmod_inf, effects = "Tave_sms", conditions = conditions,
                    re_formula = NULL) 
#compare to default priors -different!
conditional_effects(Fd_postmod_flat, effects = "Tave_sms", conditions = conditions,
                    re_formula = NULL) 

#need to extract predictions and then collapse back into ordinal categories 
library(emmeans)
epreds<-Fd_postmod_inf%>%
emmeans(~ CMIs + Tave_sms + edatope, var = "CMIs",
        at = list(edatope=unique(Fdabund$edatope),
                  CMIs = seq(-1,3, 1)),
        epred = TRUE, re_formula = NULL, allow_new_levels = TRUE) %>%emmeans::


#OTHER PACKAGE OPTIONS----
#ordinal package version
#start with intercept only model 
Fd_ord_int<-clm(newfeas_ord ~  1 , data=Fdfeas)
summary(Fd_ord_int)
exp(coef(Fd_ord_int)) #log odds scale 

Fd_ord_intr<-clmm(newfeas_ord ~  1 + (1|edatope), data=Fdfeas)
summary(Fd_ord_intr) #changes thresholds w/ group level term 
exp(coef(Fd_ord_intr)) 

Fd_ord_intr2<-clmm(newfeas_ord ~  1 + (1|edatope) + (1|BGC), data=Fdfeas) #fails to converge :/ try in brms?

Fd_ordr<-clmm(newfeas_ord ~  PC1 + (PC1|edatope), data=Fdfeas) #use climate fixed effects instead of BGC group term
summary(Fd_ordr)
ranef(Fd_ordr)
condVar(Fd_ordr)

exp(coef(Fd_ordr)) #log odds scale 


ggpred_Fd_ordr<-data.frame(ggpredict(Fd_ordr, terms = c("Tave_sm_ [-2, -1, 0, 1, 2, 3]"), bias_correction = TRUE))
colnames(ggpred_Fd_ordr)[c(1, 6)] = c("Tave_sm", "newfeas")

ggpred_Fd_ordr %>%
  mutate(newfeas = ordered(newfeas, rev(levels=rev(levels(newfeas))) %>%
  ggplot( aes(x = Tave_sm, y = predicted, fill = newfeas)) +
  geom_bar(position = "fill", stat = "identity")  + theme_minimal() + ggtitle("Probabilities of Fd Env. feas by avg.T sm")

plot(Fdfeas$CMI_~Fdfeas$newfeas_ord)
plot(Fdfeas$Tave_sm_~Fdfeas$newfeas_ord)

hist(prior_df$PC2)


#calculate ordinal response for plot data 
#feas.dat.sub<-mutate(feas.dat.sub, cover_rank=case_when(TotalA>15~1, 
#                                                        TotalA>0 & TotalA<7.5~3,
#                                                        TotalA==0 ~4, 
#                                                        TRUE~2))
#hist(feas.dat.sub$cover_rank)
#group_by(feas.dat.sub, cover_rank)%>%summarise(counts=n())#about even for 1, 2,3

#set as ordinal factor
#feas.dat.sub$cover_rank<-ordered(feas.dat.sub$cover_rank, levels = c(4, 3, 2, 1))
#str(feas.dat.sub$cover_rank)#good 


Fd_post_ord<-brm(modform_ord, Fdfeas,cores=3, chains=3, backend = "cmdstanr", threads = threading(4),  
             control = list(adapt_delta=0.99, max_treedepth = 11), 
             iter=2000, warmup = 1000, init = 0, family=  cumulative()) 
             #file= "outputs/brms/Fd_prior_mod_ord.Rdata") 
save(Fd_post_ord, file= "outputs/brms/Fd_post_mod_ord.Rdata")
summary(Fd_post_ord)

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

#run in brms
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
                                                TotalA>0 & TotalA<7.5~3,
                                                TotalA==0 ~0, 
                                                TRUE~2))
#set as ordinal factor
tree_dat$cover_rank<-ordered(tree_dat$cover_rank, levels = c(0,3, 2, 1))
str(tree_dat$cover_rank)

hist(tree_dat$TotalA)

#set priors 
priors <- c(set_prior("normal(0,4)", class = "Intercept"), 
            set_prior("normal(0, 1)", class = "sd")) #,
# set_prior("normal(0, 0.5)", class = "b"))

modform_ord<-bf(cover_rank~  
                  (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled||Species:MoistureRegime_clean) +  
                  (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled||Species:NutrientRegime_clean))

#set Moisture/Nutrient regimes to ordinal  and run as a fixed interaction with climate vars- TO DO!!
#still not sure if this is what we want...  
modform_ord<-bf(cover_rank~ (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled) * MoistureRegime_clean||Species)+
  (PPT_sm_scaled + TD_scaled + DD5_sp_scaled + Tmax_sm_scaled) * (NutrientRegime_clean||Species)

#trying for subset of spp 
tree_dat_sub<-subset(tree_dat, Species =="TSUGHET"|Species=="PSEUMEN"|Species=="PINUCON")

mod_ord<-brm(modform_ord, tree_dat_sub,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
             control = list(adapt_delta=0.99, max_treedepth = 11), 
             iter=6000, warmup = 1000, init = 0, family=  cumulative(),
             file= "outputs/brms/mod_ord3") 
summary(mod_ord3)
