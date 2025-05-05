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
#pulls in cleaned feasibility tables with BEC plot abundance and climate data (created in feas_tables.R) 
#simulates random abundance data for each feasibility class based on pre-set cutoffs
#runs gaussian & skewnormal regressions of climate on relative abundance to generate parameter
#estimates for model priors (zero and non-zero)
#runs posterior models on plot relative abundance with 1) flat (uninformed/default) and 
# 2) expert derived priors from feasibility ratings as estimated by the prior model
#compares posteriors from flat and expert informed models by BGC 

#currently this workflow is only run on Fd but can be repeated for multiple species 

#outstanding issues: 
#not currently accounting for ordinal nature of edatopic grid in the hierarchical model 
#model R2 values are low 

#load libs and data---- 
library(brms) 
#library(ordinal)
library(lme4)
#library(lmerTest)
#library(glmmTMB)
library(tidyverse)
library(remotes)
#library(NBZIMM) #zero infl gaussian
#library(bhsdtr2) #hierarchical ordinal bayesian  

#load feasibility plus abundance dataset with climate data 
load(file="data/feas_abund_clim_data.Rdata")

#create a unique variable for the 15 edatopic spaces  
feas.dat.clim$edatopex<-paste(feas.dat.clim$NutrientRegime_clean, feas.dat.clim$MoistureRegime_clean, sep="")
feas.dat.clim<- mutate(feas.dat.clim, edatope=case_when(edatopex=="C3"|edatopex=="C4"~"C34",
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
#remove Fs and 8s (Colin 12/20/24)
feas.dat.clim<-subset(feas.dat.clim, edatope!="F56"& edatope!="F3"& edatope!="DE78")
feas.dat.clim$edatope<-as.factor(feas.dat.clim$edatope)
unique(feas.dat.clim$edatope)

#set moist and nutrients as ordinal
#feas.dat.clim$MoistureRegime_clean<-ordered(feas.dat.clim$MoistureRegime_clean, levels = c(8,7,6,5,4,3,2,1,0))
#feas.dat.clim$NutrientRegime_clean<-ordered(feas.dat.clim$NutrientRegime_clean, levels = c("F", "E", "D", "C", "B", "A"))

#build prior model----
#simulate cover data with decreasing density function

simulate_data <- function(min_val, max_val, n) {
  # Create a decreasing probability density
  density_fn <- function(x) {
    1 / (x^(1))  # Decreasing power law function (adjust exponent to control steepness)
  }
  # Normalize the density to use as probabilities
  x_vals <- seq(min_val, max_val, length.out = 1000)
  prob_density <- density_fn(x_vals)
  prob_density <- prob_density / sum(prob_density)  # Normalize to sum to 1
  
  # Sample data based on the density
  sampled_vals <- sample(x_vals, n, replace = TRUE, prob = prob_density)
  
  return(sampled_vals)
}
  
#how many ratings within each class? 
knitr::kable(group_by(feas.dat.clim, newfeas_ord)%>%summarise(counts=n())) 

# Define the ranges for each class- use proposed cutoffs from Mariotte et al. 2014 (New Phyt) and 
# Simulate random data for each class within the given ranges 
#n from kable 
feas4_data <- simulate_data(0.01, 0.5, 24)
feas3_data <- simulate_data(0.51, 5, 2309)
feas2_data <- simulate_data(5.01, 12, 7364)
feas1_data <- simulate_data(12.01, 100,10996)

# Combine into one data frame
sim_data <- data.frame(
  value = c(feas1_data, feas2_data, feas3_data, feas4_data),
  newfeas = rep(c("1", "2", "3", "4"), 
              times = c(length(feas1_data),length(feas2_data), length(feas3_data), length(feas4_data))))

feas5_data<-data.frame(num = c(1 : 30), value=0, newfeas= 5)
feas5_data$num<-NULL
sim_data<- rbind(sim_data, feas5_data)
hist(sim_data$value)  

#create new df to assign simulated data 
prior_df<-feas.dat.clim
prior_df$sim_abund<-NA

# Assign simulated data values  to full dataframe based on matching class labels
for (feas in unique(prior_df$newfeas)) {
  # Find the corresponding simulated data for this class
  sim_values <- sim_data$value[sim_data$newfeas == feas]
    # Match the class rows in `existing_data` and assign the simulated values
  prior_df$sim_abund[prior_df$newfeas == feas] <- sim_values
}


#look at sim data 
hist(prior_df$sim_abund)
#compared to actual data 
hist(prior_df$TotalA)

plot(prior_df$newfeas_ord, prior_df$sim_abund) 


#run for Fd----
Fdfeas<-subset(prior_df, spp=="Fd")
str(Fdfeas)

#FIT MODEL 

#model predictors 
#PC1, PC2, PC3-> plot level climate PC axes- affect presence and abundance
#zone -> regional climate - affect presence and abundance
#subzone (bgc)-> local (subregional) climate- likely co-varies with edatope and plot level climate so currently not in model
#edatope -> site conditions- affect abundance only & interact with plot level climate 

hist(Fdfeas$sim_abund)
plot(as.factor(Fdfeas$newfeas), Fdfeas$sim_abund) 
#sqrt transform response
Fdfeas$sim_abund_sqrt<-sqrt(Fdfeas$sim_abund)
hist(Fdfeas$sim_abund_sqrt)

#brms 
#continuous response - 
Fd_priormod<-brm(bf(sim_abund_sqrt~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3||edatope)) , 
  data = Fdfeas,
  family = skew_normal(),
  chains = 2, iter = 5000, warmup = 2000, 
  control = list(adapt_delta = 0.9))
save(Fd_priormod, file= "outputs/brms/Fd_priormod_skew.Rdata")

#look at model output
summary(Fd_priormod)
ranef(Fd_priormod)
pp_check(Fd_priormod)
pp_check(Fd_priormod, type = "stat_2d")

pred <- posterior_predict(Fd_priormod)
bayesplot::ppc_dens_overlay(y = log(Fdfeas$sim_abund), 
                            yrep = log(pred[1:10,]))
bayes_R2(Fd_priormod)

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Fdfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))

conditional_effects(Fd_priormod, effects = "PC1", conditions = conditions,
                    re_formula = NULL) 

conditional_effects(Fd_priormod, effects = "PC2", conditions = conditions,
                       re_formula = NULL) 

conditional_effects(Fd_priormod, effects = "PC3", conditions = conditions,
                    re_formula = NULL) 

#posterior model
#look at data 
hist(feas.dat.clim$TotalA)
plot(feas.dat.clim$newfeas_ord, feas.dat.clim$TotalA) 
plot(feas.dat.clim$newfeas_ord, log(feas.dat.clim$TotalA+1)) #4s too high

#sqrt transform response
Fdfeas$TotalA_sqrt<-sqrt(Fdfeas$TotalA)
hist(Fdfeas$TotalA_sqrt)

#FIT MODEL first with default (flat) priors
Fd_postmod_flat<-brm(bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope)), 
                     # hu ~  CMIs + Tave_sms + zone),
                 data = Fdfeas, 
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000,
                 control = list(adapt_delta = 0.9))
                 
save(Fd_postmod_flat, file= "outputs/brms/Fd_postmod_flat.Rdata")

summary(Fd_postmod_flat)
pp_check(Fd_postmod_flat)
pred <- posterior_predict(Fd_postmod_flat)
bayesplot::ppc_dens_overlay(y = log1p(Fdabund$TotalA), 
                            yrep = log1p(pred[1:10,]))
bayes_R2(Fd_postmod_flat)
#flatpriors<-get_prior(Fd_postmod_flat)

#set informed priors from prior model
get_prior(Fd_postmod_flat)
summary(Fd_priormod)
#use mu estimates from prior model and sd estimates 
#may want to increase the width of sd estimates 
priors<- c(set_prior("normal(2.37, 0.5)", class = "Intercept"), 
          #set_prior("normal(0.45, 0.62)", class = "Intercept", dpar='hu'),
          set_prior("normal(-0.07, 0.03)", class = "b", coef = "PC1"),
          set_prior("normal(-0.14, 0.05)", class = "b", coef = "PC2"),
          set_prior("normal( 0.22, 0.12)", class = "b", coef = "PC3"),
          set_prior("normal(1.83, 0.51)", class = "b", coef = "zoneCDF"),
          set_prior("normal(1.45, 0.48)", class = "b", coef = "zoneCWH"),
          set_prior("normal(0.31, 0.48)", class = "b", coef = "zoneESSF"),
          set_prior("normal(2.16, 0.45)", class = "b", coef = "zoneICH"),
          set_prior("normal(1.4, 0.44)", class = "b", coef = "zoneIDF"),
          set_prior("normal(0.74, 0.46)", class = "b", coef = "zoneMS"),
          set_prior("normal(0.91, 0.46)", class = "b", coef = "zonePP"),
          set_prior("normal(0.79, 0.47)", class = "b", coef = "zoneSBS"),
          #set_prior("normal(-1.85, 0.14)", class = "b", coef = "Tave_sms", dpar='hu'),
          #set_prior("normal(0.83, 0.18)", class = "b", coef = "CMIs", dpar='hu'),
          set_prior("normal(0.54, 0.18)", class = "sd", group = "edatope", coef = "Intercept"), 
          set_prior("normal(0.07, 0.03)", class = "sd", group = "edatope", coef = "PC1"),   
          set_prior("normal(0.08, 0.04)", class = "sd", group = "edatope", coef = "PC2"),
          set_prior("normal(0.31, 0.12)", class = "sd", group = "edatope", coef = "PC3"),
          #set_prior("normal(2.16, 0.62)", class = "sd", group = "edatope", dpar='hu', coef = "Intercept"),  
          #set_prior("normal(0.20, 0.15)", class = "sd", group = "edatope", dpar='hu', coef = "Tave_sms"),   
          #set_prior("normal(0.32, 0.24)", class = "sd", group = "edatope", dpar='hu', coef = "CMIs"), 
          set_prior("normal(1.83,  0.03)", class="sigma"), 
          set_prior("normal(3.94,  0.24)", class="alpha")) 
          #not sure what parameters these correspond with but filling in so not run as default... 
          #set_prior("normal(0, 1)", class = "sd", lb=0),      
          #set_prior("normal(0, 1)", class = "sd", lb=0), 
          #set_prior("normal(0, 1)", class = "sd", lb=0, group = "edatope"),      
          #set_prior("normal(0, 1)", class = "sd", lb=0, group = "edatope", dpar='hu'), 
          #set_prior("normal(0, 1)", class = "b"),      
          #set_prior("normal(0, 1)", class = "b",  dpar='hu' ))      

#FIT MODEL again with expert informed priors
Fd_postmod_expert<-brm(bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope)), 
                       data = Fdfeas, prior = priors,
                       family = skew_normal(),
                       chains = 2, iter = 5000, warmup = 2000,
                       control = list(adapt_delta = 0.9),sample_prior = TRUE)
save(Fd_postmod_expert, file= "outputs/brms/Fd_postmod_expert.Rdata")
get_prior(Fd_postmod_expert)

#look at model output
summary(Fd_postmod_expert)
ranef(Fd_postmod_expert)
pp_check(Fd_postmod_expert)
pred <- posterior_predict(Fd_postmod_expert)
bayesplot::ppc_dens_overlay(y = log1p(Fdabund$TotalA), 
                            yrep = log1p(pred[1:10,]))
#not much difference in pp_check compared to post-flat prior model :/ not sure if priors are actually getting set correctly or defaults still being used?.. 12/19/24
#slightly overpredicting at lower abundances and underpredicting at higher abundances... how to update priors or sim data to address this? 

#r2
bayes_R2(Fd_postmod_expert) 
bayes_R2(Fd_postmod_flat)

#prior-posterior plots 
MODform<-bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope))
Fd_prioronly_mod<- 
  brm(MODform, cores=3, prior = priors,  data=Fdfeas, family = skew_normal(),
      sample_prior = "only") 

summary_prior<-summary(Fd_prioronly_mod)
summary_prior<-summary_prior$fixed
summary_prior$mod<-"prior"

summary_post<-summary(Fd_postmod_expert)
summary_post<-summary_post$fixed
summary_post$mod<-"posterior"

summary_all<-rbind(summary_post, summary_prior)
summary_all$param<-row.names(summary_all)
summary_all$paramx<- gsub("[^A-Za-z]+", "", summary_all$param)

#plot posterior and prior model estimates 
library(ggplot2)
ggplot(summary_all, aes(y=Estimate, x=param, fill=mod, color=mod)) +
  geom_pointrange(aes(ymin=Estimate-Est.Error, 
                      ymax=Estimate+Est.Error))+
  xlab("Model parameter")+ 
  scale_color_discrete(name="Model")+ scale_fill_discrete(name="Model")

#ok priors are not capturing posteriors... what to do now?!
plot(hypothesis(Fd_postmod_expert, "PC1 < 0"))
plot(hypothesis(Fd_postmod_expert, "PC2 > 0"))
plot(hypothesis(Fd_postmod_expert, "PC3 > 0"))
plot(hypothesis(Fd_postmod_expert, "zoneCDF> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneCWH> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneESSF< 0"))
plot(hypothesis(Fd_postmod_expert, "zoneICH> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneIDF> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneMS< 0"))
plot(hypothesis(Fd_postmod_expert, "zonePP< 0"))
plot(hypothesis(Fd_postmod_expert, "zoneSBS< 0"))

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Fdfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))
#CMI
df1<-conditional_effects(Fd_postmod_expert, effects = "PC3", conditions = conditions,
                         re_formula = NULL)
df1<-df1$PC3
df1$prior<-"expert"

#compare to default priors -different!
df2<-conditional_effects(Fd_postmod_flat, effects = "PC3", conditions = conditions,
                    re_formula = NULL) 
df2<-df2$PC3
df2$prior<-"flat"
df<-rbind(df1, df2)

ggplot(df, aes(x=PC3, y=estimate__, fill=prior, color=prior))+
         geom_smooth()+ 
         facet_wrap(~factor(edatope, levels = c("AB0", "C0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56")))+
         theme_bw()
unique(df$edatope)
#Tave
df3<-conditional_effects(Fd_postmod_expert, effects = "Tave_sms", conditions = conditions,
                    re_formula = NULL) 
df3<-df3$Tave_sms
df3$prior<-"expert"

#compare to default priors -different!
df4<-conditional_effects(Fd_postmod_flat, effects = "Tave_sms", conditions = conditions,
                    re_formula = NULL) 
df4<-df4$Tave_sms
df4$prior<-"flat"

df0<-rbind(df3, df4)

dummy<-df0[1,]
dummy2<-df0[1,]f
dummy$edatope<- "DE0"
dummy2$edatope<-"DE7"

df0<-rbind(df0, dummy, dummy2)

#back scale Temp
#mean(Fdabund$Tave_sm)
#sd(Fdabund$Tave_sm)
df0$Tave_sm<-(df0$Tave_sms+12.92)*2.27

ggplot(df0, aes(x=Tave_sm, y=estimate__, fill=prior, color=prior))+
  geom_line()+ 
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
               theme_bw() + ylab("plot rel. abundance")


#need to extract predictions and then collapse back into ordinal categories 
library(emmeans)
epreds<-Fd_postmod_expert%>%
emmeans(~ PC1 + PC2 + PC3 + edatope, var = "PC3",
        at = list(edatope=unique(Fdfeas$edatope),
                  PC3 = seq(-1,3, 1)),
        epred = TRUE, re_formula = NULL, allow_new_levels = TRUE) 



#run for Cw----
Cwfeas<-subset(prior_df, spp=="Cw")
str(Cwfeas)

#FIT MODEL 

#model predictors 
#PC1, PC2, PC3-> plot level climate PC axes- affect presence and abundance
#zone -> regional climate - affect presence and abundance
#subzone (bgc)-> local (subregional) climate- likely co-varies with edatope and plot level climate so currently not in model
#edatope -> site conditions- affect abundance only & interact with plot level climate 

hist(Cwfeas$sim_abund)
plot(as.factor(Cwfeas$newfeas), Cwfeas$sim_abund) 
#add 4s -> convert 5 to 4 for CWHvh3/103.2
Cwfeas<-mutate(Cwfeas, newfeas=if_else(newfeas==5, 4, newfeas))

#sqrt transform response
Cwfeas$sim_abund_sqrt<-sqrt(Cwfeas$sim_abund)
hist(Cwfeas$sim_abund_sqrt)

#brms 
#continuous response - 
Cw_priormod<-brm(bf(sim_abund_sqrt~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3||edatope)) , 
                 data = Cwfeas,
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000, 
                 control = list(adapt_delta = 0.9))
save(Cw_priormod, file= "outputs/brms/Cw_priormod_skew.Rdata")

#look at model output
summary(Cw_priormod)
ranef(Cw_priormod)
pp_check(Cw_priormod)
pred <- posterior_predict(Cw_priormod)
bayesplot::ppc_dens_overlay(y = log(Cwfeas$sim_abund), 
                            yrep = log(pred[1:10,]))
bayes_R2(Cw_priormod)

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Cwfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))

conditional_effects(Cw_priormod, effects = "PC1", conditions = conditions,
                    re_formula = NULL) 

conditional_effects(Cw_priormod, effects = "PC2", conditions = conditions,
                    re_formula = NULL) 

conditional_effects(Cw_priormod, effects = "PC3", conditions = conditions,
                    re_formula = NULL) 

#posterior model---- 
#look at data 
hist(feas.dat.clim$TotalA)
plot(feas.dat.clim$newfeas_ord, feas.dat.clim$TotalA) 
plot(feas.dat.clim$newfeas_ord, log(feas.dat.clim$TotalA+1)) #4s too high

#sqrt transform response
Cwfeas$TotalA_sqrt<-sqrt(Cwfeas$TotalA)
hist(Cwfeas$TotalA_sqrt)

#FIT MODEL first with default (flat) priors
Cw_postmod_flat<-brm(bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope)), 
                     # hu ~  CMIs + Tave_sms + zone),
                     data = Cwfeas, 
                     family = skew_normal(),
                     chains = 2, iter = 5000, warmup = 2000,
                     control = list(adapt_delta = 0.9))

save(Cw_postmod_flat, file= "outputs/brms/Cw_postmod_flat.Rdata")

summary(Cw_postmod_flat)
pp_check(Cw_postmod_flat)
pred <- posterior_predict(Cw_postmod_flat)
bayesplot::ppc_dens_overlay(y = log1p(Cwabund$TotalA), 
                            yrep = log1p(pred[1:10,]))
bayes_R2(Cw_postmod_flat)
#flatpriors<-get_prior(Cw_postmod_flat)

#set informed priors from prior model
get_prior(Cw_postmod_flat)
summary(Cw_priormod)
#use mu estimates from prior model and sd estimates 
#may want to increase the width of sd estimates 
priors<- c(set_prior("normal(2.37, 0.5)", class = "Intercept"), 
           #set_prior("normal(0.45, 0.62)", class = "Intercept", dpar='hu'),
           set_prior("normal(-0.07, 0.03)", class = "b", coef = "PC1"),
           set_prior("normal(-0.14, 0.05)", class = "b", coef = "PC2"),
           set_prior("normal( 0.22, 0.12)", class = "b", coef = "PC3"),
           set_prior("normal(1.83, 0.51)", class = "b", coef = "zoneCDF"),
           set_prior("normal(1.45, 0.48)", class = "b", coef = "zoneCWH"),
           set_prior("normal(0.31, 0.48)", class = "b", coef = "zoneESSF"),
           set_prior("normal(2.16, 0.45)", class = "b", coef = "zoneICH"),
           set_prior("normal(1.4, 0.44)", class = "b", coef = "zoneIDF"),
           set_prior("normal(0.74, 0.46)", class = "b", coef = "zoneMS"),
           set_prior("normal(0.91, 0.46)", class = "b", coef = "zonePP"),
           set_prior("normal(0.79, 0.47)", class = "b", coef = "zoneSBS"),
           #set_prior("normal(-1.85, 0.14)", class = "b", coef = "Tave_sms", dpar='hu'),
           #set_prior("normal(0.83, 0.18)", class = "b", coef = "CMIs", dpar='hu'),
           set_prior("normal(0.54, 0.18)", class = "sd", group = "edatope", coef = "Intercept"), 
           set_prior("normal(0.07, 0.03)", class = "sd", group = "edatope", coef = "PC1"),   
           set_prior("normal(0.08, 0.04)", class = "sd", group = "edatope", coef = "PC2"),
           set_prior("normal(0.31, 0.12)", class = "sd", group = "edatope", coef = "PC3"),
           #set_prior("normal(2.16, 0.62)", class = "sd", group = "edatope", dpar='hu', coef = "Intercept"),  
           #set_prior("normal(0.20, 0.15)", class = "sd", group = "edatope", dpar='hu', coef = "Tave_sms"),   
           #set_prior("normal(0.32, 0.24)", class = "sd", group = "edatope", dpar='hu', coef = "CMIs"), 
           set_prior("normal(1.83,  0.03)", class="sigma"), 
           set_prior("normal(3.94,  0.24)", class="alpha")) 
#not sure what parameters these correspond with but filling in so not run as default... 
#set_prior("normal(0, 1)", class = "sd", lb=0),      
#set_prior("normal(0, 1)", class = "sd", lb=0), 
#set_prior("normal(0, 1)", class = "sd", lb=0, group = "edatope"),      
#set_prior("normal(0, 1)", class = "sd", lb=0, group = "edatope", dpar='hu'), 
#set_prior("normal(0, 1)", class = "b"),      
#set_prior("normal(0, 1)", class = "b",  dpar='hu' ))      

#FIT MODEL again with expert informed priors
Cw_postmod_expert<-brm(bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope)), 
                       data = Cwfeas, prior = priors,
                       family = skew_normal(),
                       chains = 2, iter = 5000, warmup = 2000,
                       control = list(adapt_delta = 0.9),sample_prior = TRUE)
save(Cw_postmod_expert, file= "outputs/brms/Cw_postmod_expert.Rdata")
get_prior(Cw_postmod_expert)

#look at model output
summary(Cw_postmod_expert)
ranef(Cw_postmod_expert)
pp_check(Cw_postmod_expert)
pred <- posterior_predict(Cw_postmod_expert)
bayesplot::ppc_dens_overlay(y = log1p(Cwabund$TotalA), 
                            yrep = log1p(pred[1:10,]))
#not much difference in pp_check compared to post-flat prior model :/ not sure if priors are actually getting set correctly or defaults still being used?.. 12/19/24
#slightly overpredicting at lower abundances and underpredicting at higher abundances... how to update priors or sim data to address this? 

#r2
bayes_R2(Cw_postmod_expert) 
bayes_R2(Cw_postmod_flat)

#prior-posterior plots 
MODform<-bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope))
Cw_prioronly_mod<- 
  brm(MODform, cores=3, prior = priors,  data=Cwfeas, family = skew_normal(),
      sample_prior = "only") 

summary_prior<-summary(Cw_prioronly_mod)
summary_prior<-summary_prior$fixed
summary_prior$mod<-"prior"

summary_post<-summary(Cw_postmod_expert)
summary_post<-summary_post$fixed
summary_post$mod<-"posterior"

summary_all<-rbind(summary_post, summary_prior)
summary_all$param<-row.names(summary_all)
summary_all$paramx<- gsub("[^A-Za-z]+", "", summary_all$param)

#plot posterior and prior model estimates 
library(ggplot2)
ggplot(summary_all, aes(y=Estimate, x=param, fill=mod, color=mod)) +
  geom_pointrange(aes(ymin=Estimate-Est.Error, 
                      ymax=Estimate+Est.Error))+
  xlab("Model parameter")+ 
  scale_color_discrete(name="Model")+ scale_fill_discrete(name="Model")

#ok priors are not capturing posteriors... what to do now?!
plot(hypothesis(Cw_postmod_expert, "PC1 < 0"))
plot(hypothesis(Cw_postmod_expert, "PC2 > 0"))
plot(hypothesis(Cw_postmod_expert, "PC3 > 0"))
plot(hypothesis(Cw_postmod_expert, "zoneCDF> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneCWH> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneESSF< 0"))
plot(hypothesis(Cw_postmod_expert, "zoneICH> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneIDF> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneMS< 0"))
plot(hypothesis(Cw_postmod_expert, "zonePP< 0"))
plot(hypothesis(Cw_postmod_expert, "zoneSBS< 0"))

#plot conditional effects (mu+hu) 
conditions <- expand_grid(edatope = unique(Cwfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))
#CMI
df1<-conditional_effects(Cw_postmod_expert, effects = "PC3", conditions = conditions,
                         re_formula = NULL)
df1<-df1$PC3
df1$prior<-"expert"

#compare to default priors -different!
df2<-conditional_effects(Cw_postmod_flat, effects = "PC3", conditions = conditions,
                         re_formula = NULL) 
df2<-df2$PC3
df2$prior<-"flat"
df<-rbind(df1, df2)

ggplot(df, aes(x=PC3, y=estimate__, fill=prior, color=prior))+
  geom_smooth()+ 
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56")))+
  theme_bw()
unique(df$edatope)
#Tave
df3<-conditional_effects(Cw_postmod_expert, effects = "Tave_sms", conditions = conditions,
                         re_formula = NULL) 
df3<-df3$Tave_sms
df3$prior<-"expert"

#compare to default priors -different!
df4<-conditional_effects(Cw_postmod_flat, effects = "Tave_sms", conditions = conditions,
                         re_formula = NULL) 
df4<-df4$Tave_sms
df4$prior<-"flat"

df0<-rbind(df3, df4)

dummy<-df0[1,]
dummy2<-df0[1,]f
dummy$edatope<- "DE0"
dummy2$edatope<-"DE7"

df0<-rbind(df0, dummy, dummy2)

#back scale Temp
#mean(Cwabund$Tave_sm)
#sd(Cwabund$Tave_sm)
df0$Tave_sm<-(df0$Tave_sms+12.92)*2.27

ggplot(df0, aes(x=Tave_sm, y=estimate__, fill=prior, color=prior))+
  geom_line()+ 
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abundance")


#need to extract predictions and then collapse back into ordinal categories 
library(emmeans)
epreds<-Cw_postmod_expert%>%
  emmeans(~ PC1 + PC2 + PC3 + edatope, var = "PC3",
          at = list(edatope=unique(Cwfeas$edatope),
                    PC3 = seq(-1,3, 1)),
          epred = TRUE, re_formula = NULL, allow_new_levels = TRUE) 



#other package options----
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



#calculate ordinal response for plot data 
#feas.dat.clim<-mutate(feas.dat.clim, cover_rank=case_when(TotalA>15~1, 
#                                                        TotalA>0 & TotalA<7.5~3,
#                                                        TotalA==0 ~4, 
#                                                        TRUE~2))
#hist(feas.dat.clim$cover_rank)
#group_by(feas.dat.clim, cover_rank)%>%summarise(counts=n())#about even for 1, 2,3

#set as ordinal factor
#feas.dat.clim$cover_rank<-ordered(feas.dat.clim$cover_rank, levels = c(4, 3, 2, 1))
#str(feas.dat.clim$cover_rank)#good 


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
