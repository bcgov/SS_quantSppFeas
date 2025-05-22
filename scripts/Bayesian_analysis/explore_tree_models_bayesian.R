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
#simulates random abundance data for each feasibility class based on pre-set cutoffs w/overlap
#runs skewnormal regressions of plot and regional climate plus site (edatope) on simulated relative abundances 
#to generate parameter estimates for data model priors 
#runs data model  on plot relative abundance with 1) flat (uninformed/default) and 
# 2) expert derived priors from feasibility ratings as estimated by the prior model
#compares posteriors from flat and expert informed models by BGC 
#compares priors vs posteriors as cross-validation 

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

#build prior model dataset
#simulate cover data with decreasing density function
simulate_data <- function(min_val, max_val, n) {
  # Create a decreasing probability density
  density_fn <- function(x) {
    1 / (x^(1.25))  # Decreasing power law function (adjust exponent to control steepness)
  }
  # Normalize the density to use as probabilities
  x_vals <- seq(min_val, max_val, length.out = 1000)
  prob_density <- density_fn(x_vals)
  prob_density <- prob_density / sum(prob_density)  # Normalize to sum to 1
  
  # Sample data based on the density
  sampled_vals <- sample(x_vals, n, replace = TRUE, prob = prob_density)
  
  return(sampled_vals)
}

#run for Fd----
Fdfeas<-subset(feas.dat.clim, spp=="Fd")
#make sure all zones have >= 3 records 
group_by(Fdfeas, zone)%>%summarise(ct=n())

#how many ratings within each class? 
nrat<-group_by(Fdfeas, newfeas_ord)%>%summarise(counts=n())
nrat<-nrat$counts

# Define the ranges for each class- use proposed cutoffs from Mariotte et al. 2014 (New Phyt) and 
# Simulate random data for each class within the given ranges 
#feas4_data <- simulate_data(0.01, 0.5, nrat[1])
#feas3_data <- simulate_data(0.51, 2, nrat[2]) 
#feas2_data <- simulate_data(2.01, 12, nrat[3])
#feas1_data <- simulate_data(12.01, 100, nrat[4])

#same cutoffs but allowing for 20% overlap between the ranges 
feas4_data <- simulate_data(0.01, 0.8, nrat[1])
feas3_data <- simulate_data(0.51, 4, nrat[2]) 
feas2_data <- simulate_data(2.01, 30, nrat[3])
feas1_data <- simulate_data(12.01, 100, nrat[4])

#try using example area cutoffs from Table 3.3 LMH 25
#feas4_data <- simulate_data(0.01, 1, nrat[1])
#feas3_data <- simulate_data(1.01, 10, nrat[2])
#feas2_data <- simulate_data(10.01, 25, nrat[3])
#feas1_data <- simulate_data(25.01, 100, nrat[4])


# Combine into one data frame
sim_data <- data.frame(
  value = c(feas1_data, feas2_data, feas3_data, feas4_data),
  newfeas = rep(c("1", "2", "3", "4"), 
              times = c(length(feas1_data),length(feas2_data), length(feas3_data), length(feas4_data))))

#add zeroes 
#feas5_data<-data.frame(num = c(1 : 30), value=0, newfeas= 5)
#feas5_data$num<-NULL
#sim_data<- rbind(sim_data, feas5_data)
#hist(sim_data$value)  

#create new col to assign simulated data 
Fdfeas$sim_abund<-NA

# Assign simulated data values  to full dataframe based on matching class labels
for (feas in unique(Fdfeas$newfeas)) {
  # Find the corresponding simulated data for this class
  sim_values <- sim_data$value[sim_data$newfeas == feas]
    # Match the class rows in `existing_data` and assign the simulated values
 Fdfeas$sim_abund[Fdfeas$newfeas == feas] <- sim_values
}

#look at sim data 
hist(Fdfeas$sim_abund)
#look at actual data 
hist(Fdfeas$TotalA)
#by level
plot(as.factor(Fdfeas$newfeas), Fdfeas$sim_abund) 

#sqrt transform response
Fdfeas$sim_abund_sqrt<-sqrt(Fdfeas$sim_abund)
hist(Fdfeas$sim_abund_sqrt)

#fit brms models ----
#model predictors 
#zone -> regional climate  baseline zone = BG
#subzone (bgc)-> local (subregional) climate- likely co-varies with edatope and plot level climate so currently not in model
#PC1, PC2, PC3-> plot level climate (PC axes) 
#PC1~ spring/summer temp (-), PC2~ cMI (+), CMD (-), PC3~ PAS (-)
#edatope -> site conditions- SMR/SNR (interact with plot level climate)

#Prior model 
#continuous response - 
Fd_priormod<-brm(bf(sim_abund_sqrt~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3|edatope)) , 
  data = Fdfeas,
  family = skew_normal(),
  chains = 2, iter = 5000, warmup = 2000, 
  control = list(adapt_delta = 0.9))
save(Fd_priormod, file= "outputs/brms/Fd_priormod_skew.Rdata")

#look at model output
summary(Fd_priormod)
ranef(Fd_priormod)
pp_check(Fd_priormod)
pp_check(Fd_priormod, type = "stat_2d") #looks good 
bayes_R2(Fd_priormod)

#posterior model
#look at data 
hist(feas.dat.clim$TotalA)
plot(feas.dat.clim$newfeas_ord, feas.dat.clim$TotalA) 
plot(feas.dat.clim$newfeas_ord, log(feas.dat.clim$TotalA+1)) #4s too high

#sqrt transform response
hist(Fdfeas$TotalA)
Fdfeas$TotalA_sqrt<-sqrt(Fdfeas$TotalA)
hist(Fdfeas$TotalA_sqrt)

#FIT MODEL first with default (flat) priors
Fd_postmod_flat<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3|edatope)), 
                 data = Fdfeas, 
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000,
                 control = list(adapt_delta = 0.9))
                 
save(Fd_postmod_flat, file= "outputs/brms/Fd_postmod_flat.Rdata")

summary(Fd_postmod_flat)
pp_check(Fd_postmod_flat)


#set informed priors from prior model
get_prior(Fd_postmod_flat)
summary(Fd_priormod)

#use mu & error estimates from prior model  
#Intercept prior- https://discourse.mc-stan.org/t/understanding-intercept-prior-in-brms/34027/11
priors<- c(set_prior("normal(2.21, 0.56)", class = "b", coef = "Intercept"), 
          set_prior("normal(-0.11, 0.04)", class = "b", coef = "PC1"),
          set_prior("normal(-0.14, 0.06)", class = "b", coef = "PC2"),
          set_prior("normal(0.16, 0.12)", class = "b", coef = "PC3"),
          set_prior("normal(2.05, 0.57)", class = "b", coef = "zoneCDF"),
          set_prior("normal(1.82, 0.54)", class = "b", coef = "zoneCWH"),
          set_prior("normal(0.73, 0.56)", class = "b", coef = "zoneESSF"),
          set_prior("normal(2.52, 0.52)", class = "b", coef = "zoneICH"),
          set_prior("normal(1.56, 0.52)", class = "b", coef = "zoneIDF"),
          set_prior("normal(0.97, 0.54)", class = "b", coef = "zoneMS"),
          set_prior("normal(0.73, 0.55)", class = "b", coef = "zonePP"),
          set_prior("normal(0.93, 0.54)", class = "b", coef = "zoneSBS"),
          #set_prior("normal(0.77, 0.22)", class = "sd", group = "edatope"), 
          set_prior("normal(0.51, 0.19)", class = "sd", group = "edatope", coef = "Intercept"), 
          set_prior("normal(0.08, 0.03)", class = "sd", group = "edatope", coef = "PC1"),   
          set_prior("normal(0.13, 0.05)", class = "sd", group = "edatope", coef = "PC2"),
          set_prior("normal(0.27, 0.11)", class = "sd", group = "edatope", coef = "PC3"),
          set_prior("lkj(2)", class = "cor"), #set weak prior on all group correlation terms 
          set_prior("normal(2.10,  0.03)", class="sigma"), 
          set_prior("normal(3.36,  0.24)", class="alpha")) 
       
#FIT MODEL again with expert informed priors
Fd_postmod_expert<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3|edatope)), 
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
bayes_R2(Fd_postmod_expert) 

#prior-posterior plots ----
MODform<-bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 +PC2 + PC3|edatope))
Fd_prioronly_mod<- 
  brm(MODform, cores=3, prior = priors,  data=Fdfeas, family = skew_normal(),
      sample_prior = "only") 

#pop effects
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

group_effects<-ranef(Fd_postmod_expert)
group_effects<-as.data.frame(group_effects$edatope)


#plot prior-posterior checks 
plot(hypothesis(Fd_postmod_expert, "PC1 < 0"))
plot(hypothesis(Fd_postmod_expert, "PC2 > 0"))
plot(hypothesis(Fd_postmod_expert, "PC3 > 0"))
plot(hypothesis(Cw_postmod_expert, "Intercept> 0")) #captures BG
plot(hypothesis(Fd_postmod_expert, "zoneCDF> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneCWH> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneESSF< 0"))
plot(hypothesis(Fd_postmod_expert, "zoneICH> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneIDF> 0"))
plot(hypothesis(Fd_postmod_expert, "zoneMS< 0"))
plot(hypothesis(Fd_postmod_expert, "zonePP< 0"))
plot(hypothesis(Fd_postmod_expert, "zoneSBS< 0"))

#conditional effects----
conditions <- expand_grid(edatope = unique(Fdfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))
ce1<-conditional_effects(Fd_postmod_expert, effects = "PC1", conditions = conditions,
                    re_formula = NULL) 

ce2<-conditional_effects(Fd_postmod_expert, effects = "PC2", conditions = conditions,
                    re_formula = NULL) 

ce3<-conditional_effects(Fd_postmod_expert, effects = "PC3", conditions = conditions,
                    re_formula = NULL)
#plot manually
ce1<-as.data.frame(ce1$PC1)
ggplot(ce1, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC1")
ce2<-as.data.frame(ce2$PC2)
ggplot(ce2, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC2")
ce3<-as.data.frame(ce3$PC3)
ggplot(ce3, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC3")


#kfold CV----
looFd<-loo(Fd_postmod_expert, save_psis = TRUE)
kfoldFd<- kfold(Fd_postmod_expert, K = 5, save_fits = T)
save(kfoldFd, file= "outputs/brms/kfold_Fd.Rdata")
predict_Fd <- kfold_predict(kfoldFd, method = "predict")

MAE <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  mean(abs(yrep_mean - y))
}
kfold_mae <- MAE(y = predict_Fd$y, yrep = predict_Fd$yrep)
print(kfold_mae)



#intercept only models ----
Fd_priormod_int<-brm(bf(sim_abund_sqrt~ (1|zone:edatope)) , 
                     data = Fdfeas,
                     family = skew_normal(),
                     chains = 2, iter = 5000, warmup = 2000, 
                     control = list(adapt_delta = 0.9))
pp_check(Fd_priormod_int)
save(Fd_priormod_int, file= "outputs/brms/Fd_priormod_interceptonly.Rdata")
summary(Fd_priormod_int)

Fdpriors<- c(set_prior("normal(4.2, 0.14)", class= "Intercept"), 
           set_prior("normal(1.01, 0.11)", class = "sd"), 
           set_prior("normal(2.08,  0.03)", class="sigma"), 
           set_prior("normal(3.04,  0.22)", class="alpha"))

Fd_postmod_int<-brm(bf(TotalA_sqrt~ (1|zone:edatope)) , 
                    data = Fdfeas,
                    family = skew_normal(), prior =Fdpriors, 
                    chains = 2, iter = 5000, warmup = 2000, 
                    control = list(adapt_delta = 0.9))
pp_check(Fd_postmod_int)
save(Fd_postmod_int, file= "outputs/brms/Fd_postmod_interceptonly.Rdata")


#run for Cw----
Cwfeas<-subset(feas.dat.clim, spp=="Cw")
#set 5s with >0 plot abundance as 4s
Cwfeas<-mutate(Cwfeas, newfeas=ifelse(newfeas==5, 4, newfeas))

#make sure all zones have >= 3 records 
group_by(Cwfeas, zone)%>%summarise(ct=n()) 

#how many ratings within each class? 
nrat<-group_by(Cwfeas, newfeas)%>%summarise(counts=n())
nrat<-nrat$counts

# Define the ranges for each class- use proposed cutoffs from Mariotte et al. 2014 (New Phyt) and 
# Simulate random data for each class within the given ranges 
#same cutoffs but allowing for 20% overlap between the ranges 
feas4_data <- simulate_data(0.01, 0.8, nrat[4])
feas3_data <- simulate_data(0.51, 4, nrat[3]) 
feas2_data <- simulate_data(2.01, 30, nrat[2])
feas1_data <- simulate_data(12.01, 100, nrat[1])

# Combine into one data frame
sim_data <- data.frame(
  value = c(feas1_data, feas2_data, feas3_data, feas4_data),
  newfeas = rep(c("1", "2", "3", "4"), 
                times = c(length(feas1_data),length(feas2_data), length(feas3_data), length(feas4_data))))

#create new col to assign simulated data 
Cwfeas$sim_abund<-NA

# Assign simulated data values  to full dataframe based on matching class labels
for (feas in unique(Cwfeas$newfeas)) {
  # Find the corresponding simulated data for this class
  sim_values <- sim_data$value[sim_data$newfeas == feas]
  # Match the class rows in `existing_data` and assign the simulated values
  Cwfeas$sim_abund[Cwfeas$newfeas == feas] <- sim_values
}

#look at sim data 
hist(Cwfeas$sim_abund)
#compared to actual data 
hist(Cwfeas$TotalA)

plot(as.factor(Cwfeas$newfeas), Cwfeas$sim_abund) 

#sqrt transform response
Cwfeas$sim_abund_sqrt<-sqrt(Cwfeas$sim_abund)
hist(Cwfeas$sim_abund_sqrt)

#fit brms models ----
#model predictors 
#zone -> regional climate #baseline = CDF
#subzone (bgc)-> local (subregional) climate- likely co-varies with edatope and plot level climate so currently not in model
#PC1, PC2, PC3-> plot level climate (PC axes) 
#PC1~ spring/summer temp (-), PC2~ cMI (+), CMD (-), PC3~ PAS (-)
#edatope -> site conditions- SMR/SNR (interact with plot level climate)


#brms 
Cw_priormod<-brm(bf(sim_abund_sqrt~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3|edatope)) , 
                 data = Cwfeas,
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000, 
                 control = list(adapt_delta = 0.9))
save(Cw_priormod, file= "outputs/brms/Cw_priormod_skew.Rdata")

#look at model output
summary(Cw_priormod)
ranef(Cw_priormod)
pp_check(Cw_priormod)
pp_check(Cw_priormod, type = "stat_2d") #under estimating sds ~25%
bayes_R2(Cw_priormod)

#posterior model
#look at data 
hist(Cwfeas$TotalA)
#sqrt transform response
Cwfeas$TotalA_sqrt<-sqrt(Cwfeas$TotalA)
hist(Cwfeas$TotalA_sqrt)

#FIT MODEL first with default (flat) priors
Cw_postmod_flat<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 |edatope)), 
                     # hu ~  CMIs + Tave_sms + zone),
                     data = Cwfeas, 
                     family = skew_normal(),
                     chains = 2, iter = 5000, warmup = 2000,
                     control = list(adapt_delta = 0.9))

save(Cw_postmod_flat, file= "outputs/brms/Cw_postmod_flat.Rdata")

summary(Cw_postmod_flat)
pp_check(Cw_postmod_flat)
bayes_R2(Cw_postmod_flat)

#set informed priors from prior model
get_prior(Cw_postmod_flat)
summary(Cw_priormod)
#use mu estimates from prior model and sd estimates 
#may want to increase the width of sd estimates ~25%
#Intercept prior- https://discourse.mc-stan.org/t/understanding-intercept-prior-in-brms/34027/11
priors<- c(set_prior("normal(1.94, 0.40)", class = "b", coef = "Intercept"), 
           set_prior("normal(-0.09, 0.05)", class = "b", coef = "PC1"),
           set_prior("normal(0.22, 0.04)", class = "b", coef = "PC2"),
           set_prior("normal(0.10, 0.05)", class = "b", coef = "PC3"),
           set_prior("normal(1.06, 0.21)", class = "b", coef = "zoneCWH"),
           set_prior("normal(0.52, 0.32)", class = "b", coef = "zoneESSF"),
           set_prior("normal(2.48, 0.25)", class = "b", coef = "zoneICH"),
           set_prior("normal(0.09, 0.60)", class = "b", coef = "zoneIDF"),
           set_prior("normal(-1.39, 1.12)", class = "b", coef = "zoneMH"),
           #set_prior("normal(0.77, 0.22)", class = "sd", group = "edatope"), 
           set_prior("normal(0.93, 0.33)", class = "sd", group = "edatope", coef = "Intercept"), 
           set_prior("normal(0.12, 0.04)", class = "sd", group = "edatope", coef = "PC1"),   
           set_prior("normal(0.12, 0.05)", class = "sd", group = "edatope", coef = "PC2"),
           set_prior("normal(0.08, 0.06)", class = "sd", group = "edatope", coef = "PC3"),
           set_prior("lkj(2)", class = "cor"), #set weak prior on all correlation terms 
           #set_prior("normal(-0.37, 0.33)", class = "L", group = "edatope", coef = "Intercept, PC1"),
           #set_prior("normal(0.45, 0.30)", class = "L", group = "edatope", coef = "Intercept, PC2"),
           #set_prior("normal(-0.3, 0.34)", class = "L", group = "edatope", coef = "Intercept, PC3"),
           #set_prior("normal(0.07, 0.34)", class = "L", group = "edatope", coef = "PC1,PC2"),
           #set_prior("normal(0.06, 0.36)", class = "L", group = "edatope", coef = "PC2, PC3"),
           #set_prior("normal(-0.45, 0.32)", class = "L", group = "edatope", coef = "PC2, PC3"),
           set_prior("normal(2.04,  0.04)", class="sigma"), 
           set_prior("normal(3.51,  0.32)", class="alpha")) 

#FIT MODEL again with expert informed priors
Cw_postmod_expert<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 |edatope)), 
                       data = Cwfeas, prior = priors,
                       family = skew_normal(),
                       chains = 2, iter = 5000, warmup = 2000,
                       control = list(adapt_delta = 0.9),sample_prior = TRUE)
save(Cw_postmod_expert, file= "outputs/brms/Cw_postmod_expert.Rdata")

#look at model output
summary(Cw_postmod_expert)
ranef(Cw_postmod_expert)
pp_check(Cw_postmod_expert)
#r2
bayes_R2(Cw_postmod_expert) 
bayes_R2(Cw_postmod_flat)

#prior-posterior plots---- 
MODform<-bf(TotalA_sqrt ~0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 |edatope))
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

plot(hypothesis(Cw_postmod_expert, "PC1 < 0"))
plot(hypothesis(Cw_postmod_expert, "PC2 > 0"))
plot(hypothesis(Cw_postmod_expert, "PC3 > 0"))
plot(hypothesis(Cw_postmod_expert, "Intercept> 0")) #captures CDF
plot(hypothesis(Cw_postmod_expert, "zoneCWH> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneESSF< 0"))
plot(hypothesis(Cw_postmod_expert, "zoneICH> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneIDF> 0"))
plot(hypothesis(Cw_postmod_expert, "zoneMH< 0"))
unique(Cw_postmod_expert$data$zone)
#conditional effects----
conditions <- expand_grid(edatope = unique(Cwfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))
ce1<-conditional_effects(Cw_postmod_expert, effects = "PC1", conditions = conditions,
                         re_formula = NULL) 

ce2<-conditional_effects(Cw_postmod_expert, effects = "PC2", conditions = conditions,
                         re_formula = NULL) 

ce3<-conditional_effects(Cw_postmod_expert, effects = "PC3", conditions = conditions,
                         re_formula = NULL)
#plot manually
ce1<-as.data.frame(ce1$PC1)
ggplot(ce1, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC1")
ce2<-as.data.frame(ce2$PC2)
ggplot(ce2, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC2")
ce3<-as.data.frame(ce3$PC3)
ggplot(ce3, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC3")

#intercept only models ----
Cw_priormod_int<-brm(bf(sim_abund_sqrt~ (1|zone:edatope)) , 
                     data = Cwfeas,
                     family = skew_normal(),
                     chains = 2, iter = 5000, warmup = 2000, 
                     control = list(adapt_delta = 0.9))
pp_check(Cw_priormod_int)
save(Cw_priormod_int, file= "outputs/brms/Cw_priormod_interceptonly.Rdata")
summary(Cw_priormod_int)

Cwpriors<- c(set_prior("normal(3.69, 0.12)", class= "Intercept"), 
           set_prior("normal(0.6, 0.11)", class = "sd"), 
           set_prior("normal(2.26,  0.04)", class="sigma"), 
           set_prior("normal(10.07, 1.15 )", class="alpha"))

Cw_postmod_int<-brm(bf(TotalA_sqrt~ (1|zone:edatope)) , 
                    data = Cwfeas,
                    family = skew_normal(), prior =Cwpriors, 
                    chains = 2, iter = 5000, warmup = 2000, 
                    control = list(adapt_delta = 0.9))
pp_check(Cw_postmod_int)
save(Cw_postmod_int, file= "outputs/brms/Cw_postmod_interceptonly.Rdata")

#run for Pl----
Plfeas<-subset(feas.dat.clim, spp=="Pl")

#make sure all zones have >= 3 records 
group_by(Plfeas, zone)%>%summarise(ct=n()) 
Plfeas<-subset(Plfeas, zone!="CDF"& zone!="MH")

#set some 3s with low plot abundance (<0.8) as 4s
Plfeas<-mutate(Plfeas, newfeas=ifelse(newfeas==3 & TotalA<0.8, 4, newfeas))

#how many ratings within each class? 
nrat<-group_by(Plfeas, newfeas)%>%summarise(counts=n())
nrat<-nrat$counts


# Define the ranges for each class- use proposed cutoffs from Mariotte et al. 2014 (New Phyt) and 
# Simulate random data for each class within the given ranges 
#same cutoffs but allowing for 20% overlap between the ranges 
feas4_data <- simulate_data(0.01, 0.8, nrat[4])
feas3_data <- simulate_data(0.51, 4, nrat[3]) 
feas2_data <- simulate_data(2.01, 30, nrat[2])
feas1_data <- simulate_data(12.01, 100, nrat[1])

# Combine into one data frame
sim_data <- data.frame(
  value = c(feas1_data, feas2_data, feas3_data, feas4_data),
  newfeas = rep(c("1", "2", "3", "4"), 
                times = c(length(feas1_data),length(feas2_data), length(feas3_data), length(feas4_data))))

#create new col to assign simulated data 
Plfeas$sim_abund<-NA

# Assign simulated data values  to full dataframe based on matching class labels
for (feas in unique(Plfeas$newfeas)) {
  # Find the corresponding simulated data for this class
  sim_values <- sim_data$value[sim_data$newfeas == feas]
  # Match the class rows in `existing_data` and assign the simulated values
  Plfeas$sim_abund[Plfeas$newfeas == feas] <- sim_values
}

#look at sim data 
hist(Plfeas$sim_abund)
#compared to actual data 
hist(Plfeas$TotalA)

plot(as.factor(Plfeas$newfeas), Plfeas$sim_abund) 

#sqrt transform response
Plfeas$sim_abund_sqrt<-sqrt(Plfeas$sim_abund)
hist(Plfeas$sim_abund_sqrt)

#fit brms models ----
#model predictors 
#zone -> regional climate #baseline= BWBS
#subzone (bgc)-> local (subregional) climate- likely co-varies with edatope and plot level climate so currently not in model
#PC1, PC2, PC3-> plot level climate (PC axes) 
#PC1~ spring/summer temp (-), PC2~ cMI (+), CMD (-), PC3~ PAS (-)
#edatope -> site conditions- SMR/SNR (interact with plot level climate)

#prior model 
Pl_priormod<-brm(bf(sim_abund_sqrt~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3|edatope)) , 
                 data = Plfeas,
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000, 
                 control = list(adapt_delta = 0.9))
save(Pl_priormod, file= "outputs/brms/Pl_priormod_skew.Rdata")

#look at model output
summary(Pl_priormod)
ranef(Pl_priormod)
pp_check(Pl_priormod)
pp_check(Pl_priormod, type= "stat_2d")#looks good 
bayes_R2(Pl_priormod)

#posterior model
#look at data 
hist(Plfeas$TotalA)
#sqrt transform response
Plfeas$TotalA_sqrt<-sqrt(Plfeas$TotalA)
hist(Plfeas$TotalA_sqrt)

#FIT MODEL first with default (flat) priors
Pl_postmod_flat<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 |edatope)), 
                     data = Plfeas, 
                     family = skew_normal(),
                     chains = 2, iter = 5000, warmup = 2000,
                     control = list(adapt_delta = 0.9))

save(Pl_postmod_flat, file= "outputs/brms/Pl_postmod_flat.Rdata")

summary(Pl_postmod_flat)
pp_check(Pl_postmod_flat)
bayes_R2(Pl_postmod_flat)

#set informed priors from prior model
get_prior(Pl_postmod_flat)
summary(Pl_priormod)
#use mu estimates from prior model and sd estimates 
#may want to increase the width of sd estimates 
#Intercept prior- https://discourse.mc-stan.org/t/understanding-intercept-prior-in-brms/34027/11
priors<- c(set_prior("normal(4.07, 0.24)", class = "b", coef = "Intercept"), 
           set_prior("normal(-0.02, 0.04)", class = "b", coef = "PC1"),
           set_prior("normal(0.14, 0.07)", class = "b", coef = "PC2"),
           set_prior("normal(0.34, 0.12)", class = "b", coef = "PC3"),
           set_prior("normal(-1.93, 0.45)", class = "b", coef = "zoneCWH"),
           set_prior("normal(-0.71, 0.13)", class = "b", coef = "zoneESSF"),
           set_prior("normal(-0.42, 0.22)", class = "b", coef = "zoneICH"),
           set_prior("normal(-0.01, 0.19)", class = "b", coef = "zoneIDF"),
           set_prior("normal(0.58, 0.15)", class = "b", coef = "zoneMS"),
           set_prior("normal(0.43, 0.24)", class = "b", coef = "zoneSBPS"),
           set_prior("normal(1.14, 0.16)", class = "b", coef = "zoneSBS"),
           set_prior("normal(-1.28, 0.63)", class = "b", coef = "zoneSWB"),
           #set_prior("normal(0.77, 0.22)", class = "sd", group = "edatope"), 
           set_prior("normal(0.35, 0.14)", class = "sd", group = "edatope", coef = "Intercept"), 
           set_prior("normal(0.05, 0.03)", class = "sd", group = "edatope", coef = "PC1"),   
           set_prior("normal(0.16, 0.06)", class = "sd", group = "edatope", coef = "PC2"),
           set_prior("normal(0.28, 0.11)", class = "sd", group = "edatope", coef = "PC3"),
           set_prior("lkj(2)", class = "cor"), #set weak prior on all correlation terms 
           #set_prior("normal(-0.37, 0.33)", class = "L", group = "edatope", coef = "Intercept, PC1"),
           #set_prior("normal(0.45, 0.30)", class = "L", group = "edatope", coef = "Intercept, PC2"),
           #set_prior("normal(-0.3, 0.34)", class = "L", group = "edatope", coef = "Intercept, PC3"),
           #set_prior("normal(0.07, 0.34)", class = "L", group = "edatope", coef = "PC1,PC2"),
           #set_prior("normal(0.06, 0.36)", class = "L", group = "edatope", coef = "PC2, PC3"),
           #set_prior("normal(-0.45, 0.32)", class = "L", group = "edatope", coef = "PC2, PC3"),
           set_prior("normal(2.09,  0.04)", class="sigma"), 
           set_prior("normal(5.10,  0.48)", class="alpha")) #too high?

#FIT MODEL again with expert informed priors
Pl_postmod_expert<-brm(bf(TotalA_sqrt ~ 0 + Intercept + PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 |edatope)), 
                       data = Plfeas, prior = priors,
                       family = skew_normal(),
                       chains = 2, iter = 5000, warmup = 2000,
                       control = list(adapt_delta = 0.9),sample_prior = TRUE)
save(Pl_postmod_expert, file= "outputs/brms/Pl_postmod_expert.Rdata")
get_prior(Pl_postmod_expert)

#look at model output
summary(Pl_postmod_expert)
ranef(Pl_postmod_expert)
pp_check(Pl_postmod_expert)
#r2
bayes_R2(Pl_postmod_expert) 
bayes_R2(Pl_postmod_flat)

#prior-posterior plots---- 
MODform<-bf(TotalA_sqrt ~ PC1 + PC2 + PC3 + zone + (PC1 + PC2 + PC3 ||edatope))
Pl_prioronly_mod<- 
  brm(MODform, cores=3, prior = priors,  data=Plfeas, family = skew_normal(),
      sample_prior = "only") 

summary_prior<-summary(Pl_prioronly_mod)
summary_prior<-summary_prior$fixed
summary_prior$mod<-"prior"

summary_post<-summary(Pl_postmod_expert)
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

plot(hypothesis(Pl_postmod_expert, "PC1 < 0"))
plot(hypothesis(Pl_postmod_expert, "PC2 > 0"))
plot(hypothesis(Pl_postmod_expert, "PC3 > 0"))
plot(hypothesis(Pl_postmod_expert, "Intercept> 0")) #captures BWBS
plot(hypothesis(Pl_postmod_expert, "zoneCWH> 0"))
plot(hypothesis(Pl_postmod_expert, "zoneESSF< 0"))
plot(hypothesis(Pl_postmod_expert, "zoneICH> 0"))
plot(hypothesis(Pl_postmod_expert, "zoneIDF> 0"))
plot(hypothesis(Pl_postmod_expert, "zoneMS< 0"))
plot(hypothesis(Pl_postmod_expert, "zoneSBPS< 0"))
plot(hypothesis(Pl_postmod_expert, "zoneSWB< 0"))

#conditional effects----
conditions <- expand_grid(edatope = unique(Plfeas$edatope)) |> 
  mutate(cond__ = paste0(edatope))
ce1<-conditional_effects(Pl_postmod_expert, effects = "PC1", conditions = conditions,
                         re_formula = NULL) 

ce2<-conditional_effects(Pl_postmod_expert, effects = "PC2", conditions = conditions,
                         re_formula = NULL) 

ce3<-conditional_effects(Pl_postmod_expert, effects = "PC3", conditions = conditions,
                         re_formula = NULL)
#plot manually
ce1<-as.data.frame(ce1$PC1)
ggplot(ce1, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC1")
ce2<-as.data.frame(ce2$PC2)
ggplot(ce2, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC2")
ce3<-as.data.frame(ce3$PC3)
ggplot(ce3, aes(x=effect1__, y=estimate__))+ geom_line() +
  geom_ribbon(aes(ymin=lower__, ymax=upper__),  alpha=0.5, fill = "lightblue")+
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("plot rel. abund (sqrt)") + xlab("PC3")


#intercept only models ----
Pl_priormod_int<-brm(bf(sim_abund_sqrt~ (1|zone:edatope)) , 
                 data = Plfeas,
                 family = skew_normal(),
                 chains = 2, iter = 5000, warmup = 2000, 
                 control = list(adapt_delta = 0.9))
pp_check(Pl_priormod_int)
save(Pl_priormod_int, file= "outputs/brms/Pl_priormod_interceptonly.Rdata")

priors<- c(set_prior("normal(3.84, 0.11)", class = "b", coef = "Intercept"), 
           set_prior("normal(0.79, 0.09)", class = "sd"), 
           set_prior("normal(2.06,  0.04)", class="sigma"), 
           set_prior("normal(5.57,  0.46)", class="alpha"))
           
Pl_postmod_int<-brm(bf(TotalA_sqrt~ (1|zone:edatope)) , 
                     data = Plfeas,
                     family = skew_normal(), prior =priors, 
                     chains = 2, iter = 5000, warmup = 2000, 
                     control = list(adapt_delta = 0.9))
pp_check(Pl_postmod_int)
save(Pl_postmod_int, file= "outputs/brms/Pl_postmod_interceptonly.Rdata")


