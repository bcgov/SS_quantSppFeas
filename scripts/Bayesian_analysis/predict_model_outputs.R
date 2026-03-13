library(tidybayes)
library(dplyr)
library(brms)
library(bayesplot)

#load models 
load("outputs/brms/postmod_interceptonly.Rdata")
#load("outputs/brms/priormod_interceptonly.Rdata")

moddat<-postmod_int$data
#moddat2<-priormod_int$data

#predict over fitted ----
#option 1:brms posterior epred
#epreds<-posterior_epred(postmod_int, re.form = NULL)
#rm(postmod_int)
#gc()
#epreds<-as.data.frame(colMeans(epreds))
#epreds$pred_abund_cube<-epreds$`colMeans(epreds)`
#epreds$`colMeans(epreds)`<-NULL
#epreds$pred_abund<-(epreds$pred_abund_cube)^3
#save(epreds, file= "outputs/brms/postmod_epreds.Rdata")

#option 2: tidybayes epred_draws (NOT USING- maxes out RAM)
#preds <- postmod_int %>%
#  epred_draws(moddat,allow_new_levels=T,
#              re_formula = NULL) 
#preds$pred_abund<-NULL
#preds$pred_abund<-(preds$.epred)^3
#preds$TotalAB<-(preds$TotalAB_cube)^3
#save(preds, file= "outputs/brms/postmod_epred_draws.Rdata")

#load saved model preds
load("outputs/brms/postmod_epreds.Rdata")

#join back with model data 
moddat<-cbind(moddat, epreds)
moddat$TotalAB<-(moddat$TotalAB_cube)^3

#free up space 
rm(postmod_int)
rm(epreds)
gc()

#join with edatope
edat<-read.csv("data/Edatopic_v13_11.csv")
edat<-filter(edat, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))
edat<-rename(edat, edatope=Edatopic, bgc=BGC, ss_nospace=SS_NoSpace)%>%select(-Source)
edat<- mutate(edat, edatopex=case_when(edatope=="C3"|edatope=="C4"~"C34",
                                                                 edatope=="C1"|edatope=="C2"|edatope=="C0"~"C12",
                                                                 edatope=="C5"|edatope=="C6"|edatope=="C7"~"C56",
                                                                 edatope=="A3"|edatope=="A4"|edatope=="B3"|edatope=="B4"~"AB34",
                                                                 edatope=="A0"|edatope=="B0"|edatope=="A1"|edatope=="A2"|edatope=="B1"|edatope=="B2"~"AB12",
                                                                 edatope=="A7"|edatope=="B7"|edatope=="A5"|edatope=="A6"|edatope=="B5"|edatope=="B6"~"AB56",
                                                                 edatope=="D3"|edatope=="D4"|edatope=="E3"|edatope=="E4"~"DE34",
                                                                 edatope=="D0"|edatope=="E0"|edatope=="D1"|edatope=="D2"|edatope=="E1"|edatope=="E2"~"DE12",
                                                                 edatope=="D7"| edatope=="E7"|edatope=="D5"|edatope=="D6"|edatope=="E5"|edatope=="E6"~"DE56",
                                                                 TRUE~NA))
unique(edat$edatopex)
edat<-select(edat, -edatope)%>%distinct(.)

moddat<-left_join(moddat, edat)

#averages by site series 
moddat2<-group_by(moddat, spp, ss_nospace,)%>%summarise(pred_abund_ss=mean(pred_abund), TotalAB_ss= mean(TotalAB))
moddat2<-tidyr::separate(data = moddat2, col = ss_nospace, into = 'bgc', remove = F, sep = "/" )
moddat2$Zone <- gsub("[^A-Z]", "", moddat2$bgc)

#convert pred abund to esuit
moddat2<-mutate(moddat2, pred_newsuit= case_when(pred_abund_ss<1~4,
                                             pred_abund_ss>=1&pred_abund_ss<10~3,
                                             pred_abund_ss>=10&pred_abund_ss<25~2,
                                             pred_abund_ss>=25~1, TRUE~0))
                                   
#join to actual esuit
esuit<-read.csv("data/suitability_v13_25.csv")
esuit<-select(esuit, ss_nospace, spp, suitability, newsuit, mod)
esuit<-mutate(esuit, mod= case_when(mod=="HAK/SCS"|mod=="SCS/HAK"|mod=="SAS-HAK"|mod=="HK"~"SS & HK",
                                    mod=="DS & EBL"~ "DS", 
                                    mod==""~NA,
                                    mod=='cd'~"CD", TRUE ~mod))

moddat2<-left_join(moddat2, esuit)

moddat2$suit_diff<-moddat2$newsuit-moddat2$pred_newsuit
hist(moddat2$suit_diff)
#moddat2<-subset(moddat2, !is.na(suit_diff))
moddat2$X<-NULL

#visualize----
library(ggplot2)
moddat2<-tidyr::unite(moddat2, col = "spp_zone",c("spp", "Zone"), remove = F)
moddat2<-tidyr::unite(moddat2, col = "zone_spp",c("Zone", "spp"), remove = F)

ggplot(data = subset(moddat2,suit_diff<3), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~spp)+ 
  geom_hline(yintercept = 0,
    color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit')

ggplot(data = subset(moddat2,suit_diff<3), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~spp_zone)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit')

ggplot(data = subset(moddat2,suit_diff<3), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~zone_spp)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit')

#by rater 
ggplot(data = subset(moddat2,suit_diff<3), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~mod)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit')

#change_type
moddat3<-mutate(moddat2, change_type=case_when(newsuit=="1"& pred_newsuit=="2"~"E1->E2",
                                               newsuit=="1"& pred_newsuit=="3"~"E1->E3",
                                               newsuit=="1"& pred_newsuit=="4"~"E1->E4",
                                               newsuit=="2"& pred_newsuit=="1"~"E2->E1", 
                                               newsuit=="2"& pred_newsuit=="3"~"E2->E3",
                                               newsuit=="2"& pred_newsuit=="4"~"E2->E4",
                                               newsuit=="3"& pred_newsuit=="1"~"E3->E1", 
                                               newsuit=="3"& pred_newsuit=="2"~"E3->E2",
                                               newsuit=="3"& pred_newsuit=="4"~"E3->E4",
                                               newsuit=="4"& pred_newsuit=="1"~"E4->E1", 
                                               newsuit=="4"& pred_newsuit=="2"~"E4->E2",
                                               newsuit=="4"& pred_newsuit=="3"~"E4->E3",
                                               newsuit=="5"& pred_newsuit=="1"~"E5->E1", 
                                               newsuit=="5"& pred_newsuit=="2"~"E5->E2",
                                               newsuit=="5"& pred_newsuit=="3"~"E5->E3",
                                               newsuit=="5"& pred_newsuit=="4"~"E5->E4", TRUE~"no change"))
moddat3<-subset(moddat3, !is.na(suit_diff))

library(paletteer)
ggplot(moddat3, aes(x=change_type, fill=as.factor(suit_diff)))+ geom_bar()+ 
  scale_fill_paletteer_d("ggsci::default_jama") + facet_wrap(~spp)


#test set for Kiri
Sx_test<-subset(moddat3, spp=='Sx')%>%select(spp, ss_nospace, bgc, newsuit, pred_newsuit, suit_diff, change_type)
Sx_test$rule_newsuit <- round(apply(Sx_test[, c("newsuit", "pred_newsuit")], 1, median, na.rm = TRUE))
Sx_test<-mutate(Sx_test, rule_newsuit= ifelse(abs(suit_diff)>1, rule_newsuit, newsuit))
write.csv(Sx_test, "Sx_test.csv")



#predict over new -----
newdat<-tidyr::expand_grid(spp = unique(moddat$spp),
                           bgc = unique(edat$bgc), 
                           edatopex=unique(edat$edatopex),
                           StructuralStage_clean=6) #just setting this to mature forest for novel preds 

epreds_new<-posterior_epred(postmod_int, newdata = newdat, re.form = NULL, allow_new_levels= T)
epreds_new<-as.data.frame(colMeans(epreds_new))
epreds_new$pred_abund_cube<-epreds_new$`colMeans(epreds_new)`
epreds_new$`colMeans(epreds_new)`<-NULL
epreds_new$pred_abund<-(epreds_new$pred_abund_cube)^3
#save(epreds_new, file= "outputs/brms/postmod_epreds_all.Rdata")

#join back with new data 
newdat2<-cbind(newdat, epreds_new)
#join w/ ss
newdat2<-left_join(newdat2, edat)
#average by ss
newdat2<-group_by(newdat2, spp, ss_nospace)%>%summarise(pred_abund_ss=mean(pred_abund))
newdat2<-subset(newdat2, !is.na(ss_nospace)) #remove edatopic space that doesn't align to any ss 
#convert pred abund to esuit
newdat2<-mutate(newdat2, pred_newsuit= case_when(pred_abund_ss<1~4,
                                                 pred_abund_ss>=1&pred_abund_ss<10~3,
                                                 pred_abund_ss>=10&pred_abund_ss<25~2,
                                                 pred_abund_ss>=25~1, TRUE~0))
#join with suit
newdat3<-left_join(newdat2, esuit)
newdat3<-mutate(newdat3, mod=ifelse(is.na(newsuit), 'brms_pred', mod))%>%
          mutate(newsuit=ifelse(is.na(newsuit), 5, newsuit))

newdat3$suit_diff<-newdat3$newsuit-newdat3$pred_newsuit

#visualize----
ggplot(data = subset(newdat3, mod!="brms_pred"), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~mod)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit')

newdat3$Zone <- gsub("[^A-Z]", "", newdat3$ss_nospace)

ggplot(data = subset(newdat3, mod=="brms_pred"& Zone=="CWH"), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~spp)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit') + ggtitle("CWH")

ggplot(data = subset(newdat3, mod=="brms_pred"& Zone=="ESSF"), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~spp)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit') + ggtitle("ESSF")

ggplot(data = subset(newdat3, mod=="brms_pred"& Zone=="ICH"), aes(y=suit_diff))+
  geom_bar() + facet_wrap(~spp)+ 
  geom_hline(yintercept = 0,
             color='red', lty=2, alpha=0.5)+
  ylab('diff expert - model pred Esuit') + ggtitle("ICH")

#which of these have plot data? 
newdat4<-left_join(newdat3, select(moddat2, spp, ss_nospace, bgc, TotalAB_ss)) #same ~2800 without ratings in moddat 2 over fitted

#which of these have been rated at least once in same bgc- keep these novel preds? 
newdat4$bgc<-NULL
newdat4<-tidyr::separate(data = newdat4, col = ss_nospace, into = 'bgc', remove = F, sep = "/" )
newdat5<-mutate(newdat4, rated= if_else(mod!="brms_pred"|is.na(mod), 1, 0))%>% 
group_by(bgc, spp)%>%mutate(rated2=sum(rated))
newdat5<-subset(newdat5, rated2>0)
newdat5$rated<-NULL
newdat5$rated2<-NULL


#OLD----
#pull out average abundances by site series-RAW DATA
avgs<-select(moddat, bgc, ss_nospace, Species, spp, newsuit_ord, mean_abund_ss, nplots_ss)%>% distinct(.)
avgs$edatope<-NULL
avgs<-distinct(avgs)

#pull out average abundances by site series-SIM DATA
avgs3<-select(moddat, bgc, ss_nospace, Species, spp, sim_abund)%>%
  group_by(spp, ss_nospace)%>%
  summarise(sim_mean_abund_ss=mean(sim_abund, na.rm = T), sim_sd_abund_ss=sd(sim_abund, na.rm = T))
avgs3$edatope<-NULL
avgs3<-distinct(avgs3)

updated_avgs<-left_join(avgs, avgs2)
updated_avgs<-left_join(updated_avgs, avgs3)
updated_avgs$newsuit<-as.numeric(as.character(updated_avgs$newsuit_ord))
updated_avgs$diff<-updated_avgs$newsuit-updated_avgs$pred_newsuit  
write.csv(updated_avgs, "outputs/expert_modelpred_ratings.csv")
ggplot(data = updated_avgs, aes(y=diff))+ geom_bar() + facet_wrap(~spp)+ geom_hline(yintercept=0, color='red', lty=2, alpha=0.5)+ylab('diff expert - model pred Esuit')

#coloredat#color by fitted or predicted 
#check<-as.data.frame(ranef(Fd_postmod_int))
#check$zone_edatope<-row.names(check)
#check<-separate(check, zone_edatope, into = c("zone", "edatope"), sep = "_")%>%subset(zone=="ICH")
#check<-unique(check$edatope)
#preds<-mutate(preds, pred_type= if_else(edatope %in% check, "fitted", "predicted"))

#plot predictions- back transform to %

#expert model
ggplot(subset(preds_all,bgc=="CDFmm"), aes(x=(.epred)^3, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(edatope))+
    #levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Fd in CDFmm") #+ theme(legend.position = 'none') 

ggplot(subset(preds_all,bgc=="BGxw1" & !is.na(ss_nospace)), aes(x=(.epred)^3, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(ss_nospace), scales="free_x")+
  #levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Fd in BGxw1") #+ theme(legend.position = 'none') 

ggplot(subset(preds_all,bgc=="BGxw1"), aes(x=(.epred)^3, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(edatope))+
  #levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Fd in BGxw1") #+ theme(legend.position = 'none') 



 
Pl_ESSF_newdata = expand_grid(edatope = unique(Plfeas$edatope),
                      zone = "ESSF", 
                      PC1=quantile(Fdfeas$PC1), 
                      PC2=quantile(Fdfeas$PC2), 
                      PC3=quantile(Fdfeas$PC3))
Pl_ESSF_preds <- Pl_postmod_expert %>%
  epred_draws(Pl_ESSF_newdata,allow_new_levels=T,
              re_formula = NULL) 


ggplot(Pl_ESSF_preds, aes(x = PC3, y = .epred)) + stat_lineribbon(alpha = 0.5) +
  ggtitle("ESSF")+  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + ylab("pred rel. abundance")

PosteriorsPl <- as_tibble(posterior_samples(Pl_postmod_expert))

 

#Pl wet zones in ESSF ----
#Pine on mesic and wetter sites in the ESSF #265
#Question from Erica Liles: you pretty much don’t find pine in the mature forest on these sites and they aren’t
#in the veg tables in BEC guides for Skeena, yet model and experts have been ranking them a 2 sometimes. 
#how to reconcile?

#use intercept only model 
Pl_ESSF_newdata = expand_grid(zone = "ESSF", edatope= unique(Plfeas$edatope))

#add missing edatopes - DE7 wet/rich, DE0 dry/rich
edatope_miss<-data.frame(edatope="DE0", zone="ESSF")
edatope_miss<-rbind(edatope_miss, data.frame(edatope="DE7", zone="ESSF"))
Pl_ESSF_newdata<-rbind(Pl_ESSF_newdata, edatope_miss)

#predict from intercept only model
Pl_ESSF_preds <- Pl_postmod_int %>%
  epred_draws(Pl_ESSF_newdata,allow_new_levels=T,
              re_formula = NULL) 
Pl_ESSF_preds$mod<-"expert informed"

Pl_ESSF_preds2 <- Pl_priormod_int %>%
  epred_draws(Pl_ESSF_newdata,allow_new_levels=T,
              re_formula = NULL) 
Pl_ESSF_preds2$mod<-"prior"

Pl_ESSF_preds<-rbind(Pl_ESSF_preds, Pl_ESSF_preds2)

#plot predictions- back transform to %
ggplot(Pl_ESSF_preds, aes(x=(.epred)^2, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Pl in ESSF ") + theme(legend.position = 'none')


#Cw disagreement in CWH
Cw_CWH_newdata = expand_grid(zone = "CWH", edatope= unique(Cwfeas$edatope))

#add missing edatopes - DE7 wet/rich, DE0 dry/rich
edatope_miss<-data.frame(edatope="DE0", zone="CWH")
edatope_miss<-rbind(edatope_miss, data.frame(edatope="DE7", zone="CWH"))
edatope_miss<-rbind(edatope_miss, data.frame(edatope="C0", zone="CWH"))

Cw_CWH_newdata<-rbind(Cw_CWH_newdata, edatope_miss)

#predict from intercept only model
Cw_CWH_preds <- Cw_postmod_int %>%
  epred_draws(Cw_CWH_newdata,allow_new_levels=T,
              re_formula = NULL) 
Cw_CWH_preds$mod<-"expert informed"

Cw_CWH_preds2 <- Cw_priormod_int %>%
  epred_draws(Cw_CWH_newdata,allow_new_levels=T,
              re_formula = NULL) 
Cw_CWH_preds2$mod<-"prior"

Cw_CWH_preds<-rbind(Cw_CWH_preds, Cw_CWH_preds2)

#Cwot predictions- back transform to %
ggplot(Cw_CWH_preds, aes(x=(.epred)^2, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(edatope, levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Cw in CWH ") #+ theme(legend.position = 'none')


#create bgc x edat table for preds
#for novel preds- use expand grid-LATER
#newdata = expand_grid(edatope = unique(moddat2$edatope),
#                     bgc = unique(moddat2$bgc))
#newdata<-left_join(newdata, edat)
#newdata<-(newdata, !is.na(ss_nospace))%>%select(-ss_nospace)%>%distinct(.) #remove things that aren't in edatopic table 

#try preds for one spp 
#newdata$spp<-"Cw"

#add missing edatope 
#edatope_miss<-data.frame(edatope="E0", bgc = unique(Fdfeas$bgc))
#newdata<-rbind(newdata, edatope_miss)


#preds$mod<-'data + expert prior'
#preds2 <- priormod_int %>%
#  epred_draws(moddat2,allow_new_levels=T,
#              re_formula = NULL) 
#preds2$mod<-'expert prior'

#preds_all<-rbind(preds, preds2)
