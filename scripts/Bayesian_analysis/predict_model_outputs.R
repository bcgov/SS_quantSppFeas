library(marginaleffects)
library(tidybayes)
library(dplyr)

# Fd ---
#use intercept only model 
newdata = expand_grid(edatope = unique(Fdfeas$edatope),
                      bgc = unique(Fdfeas$bgc))
#add missing edatope 
edatope_miss<-data.frame(edatope="E0", bgc = unique(Fdfeas$bgc))
newdata<-rbind(newdata, edatope_miss)

preds <- Fd_postmod_int %>%
  epred_draws(newdata,allow_new_levels=T,
              re_formula = NULL) 
preds$mod<-'data + expert prior'
preds2 <- Fd_priormod_int %>%
  epred_draws(newdata,allow_new_levels=T,
              re_formula = NULL) 
preds2$mod<-'expert prior'

preds_all<-rbind(preds, preds2)

#join with edatopic table to convert to site series 
edat<-read.csv("data/Edatopic_v13_11.csv")
edat<-rename(edat, edatope=Edatopic, bgc=BGC, ss_nospace=SS_NoSpace)%>%select(-Source)
preds_all<-left_join(preds_all, edat)
preds_all$pred_abund<-preds_all$.epred^2

#calculate average abundances by site series-RAW DATA
avgs<-group_by(feas.dat.sub, bgc, ss_nospace, Species, spp, newsuit_ord)%>%
summarise(mean_abund_ss=mean(TotalAB, na.rm = T), sd_abund_ss=sd(TotalAB, na.rm = T), nplots_ss=n())
avgs<-mutate(avgs, sd_abund_ss =replace(sd_abund_ss, is.na(sd_abund_ss), 0))
avgsFd<-subset(avgs, spp=="Fd")

#calculate average abundances by site series-MODEL PREDS
avgs2<-group_by(preds, ss_nospace)%>%
  summarise(pred_mean_abund_ss=mean(pred_abund, na.rm = T), pred_sd_abund_ss=sd(pred_abund, na.rm = T))
updated_avgs<-left_join(avgsFd, avgs2)
  
#coloredat#color by fitted or predicted 
#check<-as.data.frame(ranef(Fd_postmod_int))
#check$zone_edatope<-row.names(check)
#check<-separate(check, zone_edatope, into = c("zone", "edatope"), sep = "_")%>%subset(zone=="ICH")
#check<-unique(check$edatope)
#preds<-mutate(preds, pred_type= if_else(edatope %in% check, "fitted", "predicted"))

#plot predictions- back transform to %

#expert model
ggplot(subset(preds_all,bgc=="CDFmm"), aes(x=(.epred)^2, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(edatope))+
    #levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Fd in CDFmm") #+ theme(legend.position = 'none') 

ggplot(subset(preds_all,bgc=="BGxw1" & !is.na(ss_nospace)), aes(x=(.epred)^2, color=mod)) + geom_density(alpha=0.6) +
  facet_wrap(~factor(ss_nospace), scales="free_x")+
  #levels = c("AB0", "C0", "DE0", "AB12", "C12", "DE12", "AB34", "C34", "DE34", "AB56", "C56", "DE56", "AB7", "C7", "DE7")), ncol=3)+
  theme_bw() + xlab("pred rel. abundance (%)") + ggtitle("Fd in BGxw1") #+ theme(legend.position = 'none') 

ggplot(subset(preds_all,bgc=="BGxw1"), aes(x=(.epred)^2, color=mod)) + geom_density(alpha=0.6) +
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



