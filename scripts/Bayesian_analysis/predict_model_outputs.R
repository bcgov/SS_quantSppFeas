library(marginaleffects)
library(tidybayes)
library(dplyr)

#use intercept only model 
#create bgc x edat table for preds

#predict over data review=T
newdata<-select(moddat, bgc, edatope, spp)%>%distinct(.) 

#for novel preds- use expand grid-LATER
#newdata = expand_grid(edatope = unique(moddat2$edatope),
 #                     bgc = unique(moddat2$bgc))

edat<-read.csv("data/Edatopic_v13_11.csv")
edat<-rename(edat, edatope=Edatopic, bgc=BGC, ss_nospace=SS_NoSpace)%>%select(-Source)
#newdata<-left_join(newdata, edat)
#newdata<-(newdata, !is.na(ss_nospace))%>%select(-ss_nospace)%>%distinct(.) #remove things that aren't in edatopic table 

#try preds for one spp 
#newdata$spp<-"Cw"

#add missing edatope 
#edatope_miss<-data.frame(edatope="E0", bgc = unique(Fdfeas$bgc))
#newdata<-rbind(newdata, edatope_miss)

preds <- postmod_int %>%
  epred_draws(newdata,allow_new_levels=T,
              re_formula = NULL) 
preds$mod<-'data + expert prior'
preds2 <- priormod_int %>%
  epred_draws(newdata,allow_new_levels=T,
              re_formula = NULL) 
preds2$mod<-'expert prior'

#preds_all<-rbind(preds, preds2)


#join with edatopic table to convert to site series 
preds<-left_join(preds, edat)
preds$pred_abund<-exp(preds$.epred)
hist(preds$pred_abund)

preds2<-left_join(preds2, edat)
preds2$pred_abund<-exp(preds2$.epred)
hist(preds2$pred_abund)

#pull out average abundances by site series-RAW DATA
avgs<-select(moddat, bgc, ss_nospace, Species, spp, newsuit_ord, mean_abund_ss,sd_abund_ss, nplots_ss)%>% distinct(.)
avgs$edatope<-NULL
avgs<-distinct(avgs)

#calculate average abundances by site series-MODEL PREDS
avgs2<-group_by(preds, ss_nospace)%>%
  summarise(pred_mean_abund_ss=mean(pred_abund, na.rm = T), pred_sd_abund_ss=sd(pred_abund, na.rm = T))
#avgs2$pred_mean_abund_ss<-round(avgs2$pred_mean_abund_ss) #round up to next % point
avgs2<-mutate(avgs2, pred_newsuit= case_when(pred_mean_abund_ss<1~4,
                                             pred_mean_abund_ss>=1&pred_mean_abund_ss<10~3,
                                             pred_mean_abund_ss>=10&pred_mean_abund_ss<25~2,
                                             pred_mean_abund_ss>=25~1, TRUE~0))
                                        
updated_avgs<-left_join(avgs, avgs2)
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



