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

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

#brms----
#QA QC model----
load(file = 'outputs/brms/mod_ord.rds')
moddat<-mod_ord$data
mod.all<-mod_ord

#summary----
summary(mod.all) #rhats and bulk ess all good 
bayes_R2(mod_ord) 


#pull out Unique Species x SMR/SNR slopes----
slopes<-coef(object = mod.all)
SMRslopes<-as.data.frame(slopes$`Species:MoistureRegime_clean`)
SNRslopes<-as.data.frame(slopes$`Species:NutrientRegime_clean`)

SMRslopes$SMR<-row.names(SMRslopes)
SMRslopes<-separate(SMRslopes,col =  SMR, into = c("Species", "SMR"))
SNRslopes$SNR<-row.names(SNRslopes)
SNRslopes<-separate(SNRslopes,col =  SNR, into = c("Species", "SNR"))

#plot random slopes
SMRslopes<-subset(SMRslopes, SMR!="NA")
ggplot(data=SMRslopes, aes(y = as.factor(SMR) , x=Estimate.PPT_sm_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.PPT_sm_scaled, xmax=Q97.5.PPT_sm_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
 ylab("Soil Moisture Regime")

ggplot(data=SMRslopes, aes(y = as.factor(SMR) , x=Estimate.Tmax_sm_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.Tmax_sm_scaled, xmax=Q97.5.Tmax_sm_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
   ylab("Soil Moisture Regime")

ggplot(data=SMRslopes, aes(y = as.factor(SMR) , x=Estimate.DD5_sp_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.DD5_sp_scaled, xmax=Q97.5.DD5_sp_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  #theme(axis.text.y=element_blank(),  #remove y axis labels
  #           axis.ticks.y=element_blank()) +
  ylab("Soil Moisture Regime")

ggplot(data=SMRslopes, aes(y = as.factor(SMR) , x=Estimate.TD_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.TD_scaled, xmax=Q97.5.TD_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  #theme(axis.text.y=element_blank(),  #remove y axis labels
  #           axis.ticks.y=element_blank()) +
  ylab("Soil Moisture Regime")

SNRslopes<-subset(SNRslopes, SNR!="NA")
ggplot(data=SNRslopes, aes(y = as.factor(SNR) , x=Estimate.PPT_sm_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.PPT_sm_scaled, xmax=Q97.5.PPT_sm_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  ylab("Soil Nutrient Regime")

ggplot(data=SNRslopes, aes(y = as.factor(SNR) , x=Estimate.Tmax_sm_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.Tmax_sm_scaled, xmax=Q97.5.Tmax_sm_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  ylab("Soil Nutrient Regime")

ggplot(data=SNRslopes, aes(y = as.factor(SNR) , x=Estimate.DD5_sp_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.DD5_sp_scaled, xmax=Q97.5.DD5_sp_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  #theme(axis.text.y=element_blank(),  #remove y axis labels
  #           axis.ticks.y=element_blank()) +
  ylab("Soil Nutrient Regime")

ggplot(data=SNRslopes, aes(y = as.factor(SNR) , x=Estimate.TD_scaled)) + 
  geom_pointrange(aes(xmin=Q2.5.TD_scaled, xmax=Q97.5.TD_scaled), size=0.01, alpha=0.5)+
  facet_wrap(~Species)+theme_bw()+
  geom_vline(xintercept = 0, color='black', lty=2)+
  #theme(axis.text.y=element_blank(),  #remove y axis labels
  #           axis.ticks.y=element_blank()) +
  ylab("Soil Nutrient Regime")


#which species x snr/smr are missing? 
#ABIELAS, BETUPAP, LARIOCC, PICENE, PICEGLA, PICEMAR, PICESIT, PINUCON, PINUPON, POPUTRE, PSEUMEN - 8
#LARIOCC, PINUPON, POPUTRE - 7
#PINUPON - 6 
#LARIOCC, PICEGLA, PICEGLA -0

#ABIEAMA, BETUPAP, CALLNOO, LARIOCC, PICEGLA, PICEMAR, PINUPON, POPUTRE, PSEUMEN, THUJPLI, TSUGMER - F
# PICEMAR - E

#run ppchecks ----
ppc<-pp_check(mod_ord) #looks good
#save(ppc, file='outputs/brms/ppcheck_modallspp.Rdata')

#plot for predictive accuracy----
#load(file = "outputs/brms/ppcheck_modallspp.Rdata")
ppcdata<-ppc$data
obsdata<-subset(ppcdata, is_y==TRUE)
preddata<-subset(ppcdata, is_y==FALSE) #10 predictions for every 1 observation 

obsdatax<-select(obsdata, y_id, value)%>%rename(value0=value)
  
alldata<-left_join(obsdatax, preddata)

#create confusion matrix
alldata$value0<-as.factor(alldata$value0)
alldata$value<-as.factor(alldata$value)

CM<-caret::confusionMatrix(data=alldata$value, reference=alldata$value0) # Accuracy : 0.3973 
cmtab<-as.data.frame(CM$table)
ggplot(data=cmtab, aes(y = Prediction, x=Reference, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
   theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#make plot nicer
#with stat lineribbon 
#png(file = "figures/ppc_plot3.png", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=alldata, aes(y = value, x=value0)) + 
  geom_point()+ 
  #stat_lineribbon()+
  ylab("Predicted") + xlab("Observed") +
  scale_fill_brewer() + theme_bw()
#dev.off()

#as gradient - no points 
#png(file = "figures/ppc_plot.png", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=alldata, aes(y = value, x=value0, fill = after_stat(.width))) + 
  #geom_point()+ 
  stat_lineribbon(.width = ppoints(50))+
  ylab("Predicted") + xlab("Observed") +
  scale_fill_distiller() + theme_bw()
#dev.off()

#png(file = "figures/ppc_plot4.png", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=alldata, aes(y = value, x=value0)) + 
  #geom_point()+ 
  stat_summary(
    geom = "smooth",
    fun.data = mean_cl_normal,
    fun.args = list(conf.int = 0.95),
    #group = 1,
    alpha = .5,
    color = "black",
    se = TRUE) +
  ylab("Predicted") + xlab("Observed") +
  #scale_fill_distiller() + 
  theme_bw()
#dev.off()



#plot by groups of interest??
#raw data from model for plotting
load(file = 'outputs/Aug2023/Q1_rawdata_full.Rdata')
MODDAT$y_id<-as.integer(row.names(MODDAT))

alldata<-left_join(alldata, select(MODDAT, y_id, interaction_present, CLIMATE1, TAXA1, 
                                    RESOLVED.TAXA1, UNIQUE.PAIR.ID))
png(file = "figures/ppc_plot6.png", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=alldatax, aes(y = value, x=value0)) + 
  geom_point(alpha=0.3)+ 
  geom_smooth(method="lm")+
  ylab("Predicted") + xlab("Observed") +
  facet_wrap(~CLIMATE1 + TAXA1, scales="free")+ 
  #scale_fill_brewer() + 
  theme_bw()
dev.off()

#run correlation tests 
ct<-cor.test(alldata$value, alldata$value0) #POOR CORRELATION ~0.15
brms::bayes_R2(mod.all) #12%

#predictive accuracy----

#differences on raw values then calculate mean, sd on diffs
#alldatay<- mutate(alldata, diff=abs(value0-value), accuracy=(abs(value0-diff)/value0))%>%group_by(y_id)%>%
#  mutate(mean=mean(accuracy), sd=sd(accuracy))

#pull out means and sds for each observed value 
#for every observation (y_id) -492k values 
pred_estimates_yid<-group_by(alldata, y_id) %>%mutate(mean_pred=mean(value), sd_pred=sd(value))%>%
  select(value0, mean_pred, sd_pred)%>%distinct(.)%>% mutate(diff=abs(value0-mean_pred))
#combine with metatdata
pred_estimates_yid<-left_join(pred_estimates_yid, select(MODDAT, y_id, interaction_present, CLIMATE1, TAXA1, 
                                   RESOLVED.TAXA1, UNIQUE.PAIR.ID))
#for every genus pair (unique pair id) -22k values 
pred_estimates_pairid<-mutate(alldata, diff=abs(value0-value))%>%
  group_by(UNIQUE.PAIR.ID) %>%
  mutate(mean_obs=mean(value0),mean_pred=mean(value), sd_pred=sd(value), mean_diff=mean(diff), sd_diff=sd(diff))%>%
  select(mean_obs, mean_pred, sd_pred, mean_diff, sd_diff)%>%
  distinct(.)
#combine with metatdata 
pred_estimates_pairid<-left_join(pred_estimates_pairid, select(MODDAT, interaction_present, CLIMATE1, TAXA1, 
                                                         RESOLVED.TAXA1, UNIQUE.PAIR.ID))%>%
  distinct(.)

#or unique values across genera - 6025 values 
pred_estimates<-group_by(alldata, value0) %>%mutate(mean_pred=mean(value), sd_pred=sd(value))%>%
select(value0, mean_pred, sd_pred)%>%distinct(.)%>%mutate(diff=abs(value0-mean_pred))

save(alldata, pred_estimates, pred_estimates_yid, pred_estimates_pairid, file='outputs/Sep2023/Q1_ppc_data.Rdata')

load(file = "outputs/Sep2023/Q1_ppc_data.Rdata")

#plot with 1-1 line
pdf(file = "figures/ppc_plot5.pdf", width = 8, height = 8)
ggplot(data = pred_estimates, aes(x=value0, y=mean_pred))+
  geom_abline(slope=1, intercept=0, color="red")+
  geom_pointrange(aes(ymin=mean_pred-sd_pred, ymax=mean_pred+sd_pred), alpha=0.2)+  
  ylab("Predicted") + xlab("Observed")  
  #scale_colour_gradient2(
  #  low = "red",
  #  mid = "white",
  #  high = "red",
  #  midpoint = 0,
  #  space = "Lab",
  #  na.value = "grey50",
  #  guide = "colourbar",
  #  aesthetics = "colour")+
  # theme_bw() #+theme(aes(legend.position="none"))+
  dev.off()
  

#run correlation tests on pred means
ct2<-cor.test(pred_estimates$mean_pred, pred_estimates$value0)
  
#plot of corr test on means 
library(RColorBrewer)
brewer.pal(n=8,"Set2")#get some hex codes

pdf(file = "figures/ppc_plot0.pdf", width = 8, height = 6)
ggplot(data = pred_estimates, aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="darkblue", lty=2)+
  geom_point(alpha=1)+ 
  ylab("Mean predicted log change in abundance") + xlab("Observed log change in abundance")  +
  theme_bw()+
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#ff0000",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", name="Residuals")+
  annotate("text", label="R2 = 0.38
              Pearson's R= 0.62
              95% CI (0.608, 0.639)
              t = 61.869, df = 6023
              p-value < 2.2e-16", x=5.2, y=-1.2, size=3, hjust=1)  
#+theme(aes(legend.position="none"))+

dev.off()
#look at distribution of residuals 
hist(pred_estimates$diff)
check<-subset(pred_estimates, diff<1) #5438/6025 ~90%
check2<-subset(pred_estimates, diff<0.5)#4296/6025 ~71% 
check2.1<-subset(pred_estimates, diff<0.25)#2870/6025 ~47% 
check2.2<-subset(pred_estimates, diff<0.1)#1350/6025 ~22% 

#only plot those with diff < 0.5 (71% of values) 
ggplot(subset(pred_estimates,diff<0.5), aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="darkblue", lty=2)+
  geom_point(alpha=1)+
  ylab(" Mean predicted log change in abundance") + xlab("Observed log change in abundance")  +
  theme_bw()+
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#FC8D62",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", name="Residuals")

#ok looks like predictive accuracy is only very good 
#between observed log change -1, 1

hist(check2$value0)#majority of these have changes <+/- 0.5 (logged)
hist(check2.1$value0)#majority of these have changes <+/- 0.5 (logged)

check3<-subset(check2.1, check2.1$value0<0.5)#2841/2870 ~99%
check4<-subset(check2, check2$value0<0.25) #2732/2870 ~95%

#make pretty 
hist<-ggplot(check2, aes(x=value0)) + 
  geom_histogram(aes(y=..count..), colour="black", fill="#8da0cb")+
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#8DA0CB") +
  theme_bw()+
  xlab("observed log change in abundance")+
  ylab("number of observations")+ xlim(-1,1)


####Put all together into one multipanel fig 
#Fig 2

a<-ggplot(data = pred_estimates, aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="darkblue", lty=2)+
  geom_point(alpha=1)+ 
  ylab("Predicted log change in abundance") + xlab("Observed log change in abundance")  +
  theme_bw()+
  theme(axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12))+
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#ff0000",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", name="Residuals")+
  annotate("text", label="R2 = 0.38
              Pearson's R= 0.62
              95% CI (0.608, 0.639)
              t = 61.869, df = 6023
              p-value < 2.2e-16", x=5.2, y=-1.2, size=3, hjust=1)  

b<-ggplot(subset(pred_estimates,diff<0.25), aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="black", lty=2)+
  geom_point(alpha=1)+
  ylab("Predicted log change in abundance") + 
  xlab("Observed log change in abundance")  +
  theme_bw()+ 
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#Ff0000",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", name="Residuals")+
    theme(axis.title = element_text(size = 10), 
          legend.title = element_text(size = 10))

c<-ggplot(check2.1, aes(x=value0)) + 
  geom_histogram(aes(y=..count..), colour="black", fill="#8da0cb")+
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#8DA0CB") +
  theme_bw()+
  xlab("Observed log change in abundance")+
  ylab("Number of observations")+ xlim(-1,1)+
  theme(axis.title = element_text(size = 10))


library(gridExtra)
pdf(file = "figures/Figure2.pdf", width = 12, height = 8)
grid.arrange(arrangeGrob(a, ncol=1, nrow=1), heights=c(6,2.5), widths=c(3,2), 
             arrangeGrob(b,c, ncol=1, nrow=2))
dev.off()



#prior posterior plots ----
#run prior only model 
mod.all3.prior <- 
  brm(modform3, cores=3, prior = priors,  data=tree_dat, 
      sample_prior = "only")
#      file="Revision 1 ecography/output/meta model/SEspecies_prior.rmd")

summary_prior<-summary(mod.all3.prior)
summary_prior<-summary_prior$fixed
#summary_prior<-summary_prior$spec_pars
summary_prior$mod<-"prior"

summary_post<-summary(mod.all3)
summary_post<-summary_post$fixed
#summary_post<-summary_post$spec_pars
summary_post$mod<-"posterior"

summary_all<-rbind(summary_post, summary_prior)
summary_all$param<-row.names(summary_all)
summary_all$paramx<- gsub("[^A-Za-z]+", "", summary_all$param)

#plot posterior and prior model estimates 
library(ggplot2)
ggplot(summary_all, aes(y=Estimate, x=paramx, fill=mod, color=mod)) +
  geom_pointrange(aes(ymin=Estimate-Est.Error, 
                      ymax=Estimate+Est.Error))+ 
  geom_jitter()+
  xlab("Model parameter")+ 
  scale_color_discrete(name="Model")+ scale_fill_discrete(name="Model")

# Check the normal distribution of random effects----
qqnorm(slopes$Estimate.Prop.Change.Gn2, main = "Normal Q-Q plot of random slopes",  bty="n")
qqline(slopes$Estimate.Prop.Change.Gn2, col = 2, lwd = 3, lty = 3) # a bit off on the tails 


#Kfold cross validation----
#https://avehtari.github.io/modelselection/rats_kcv.html#5_K-fold_cross-validation
#use kfold split grouped to preserves group structure within folds but assigns groups randomly to folds 
grps<-loo::kfold_split_grouped(K = 10, x = moddat$UNIQUE.PAIR.ID)
print(grps)
tab<-table(moddat$UNIQUE.PAIR.ID, grps)
colSums(tab)# look pretty even ~49k each 

options(future.globals.maxSize = 8000 * 1024^2)# give a lot of space so kf doesn't crash 


kf<-kfold(mod, folds = grps, save_fits = T ) #running very slowly


#loo-CV subsample- 
#Keep getting ERROR: can't allocate vector this size for both of these
#https://mc-stan.org/loo/articles/loo2-large-data.html
#try this? 
loo_ss_1 <-
  brms::loo_subsample(mod)
print(loo_ss_1)

r2<-loo_R2(mod)






#plot as histogram 
library(RColorBrewer)
brewer.pal(n=5,"Set2")#get some hex codes 

hist<-ggplot(slopes, aes(x=Estimate.Prop.Change.Gn2)) + 
  geom_histogram(aes(y=..count..), colour="black", fill="#8DA0CB")+
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  #geom_density(alpha=.2, fill="#8DA0CB") +
   theme_bw()+
  xlab("Strength of association")+
  ylab("Number of pairs")+
  geom_vline(aes(xintercept=0),
             color="black", linetype="dashed", size=1)

#combine 
gridExtra::grid.arrange(hist, bars, nrow=1)

#add posteriors on top 
#these are equivalent 
postmod<-as_tibble(fixef(mod, summary=F))
postmod2<-as_tibble(posterior_samples(mod, pars = "b_Prop.Change.Gn2"))

transparent_theme <- theme( axis.title.y = element_blank(), axis.title.x = element_blank(),
                             axis.text.x = element_blank(), 
                             axis.text.y = element_blank(), 
                             axis.ticks.x = element_blank(), 
                             axis.ticks.y = element_blank(), 
                             panel.grid = element_blank(), 
                             axis.line = element_blank(), 
                             panel.background = element_blank(), 
                             plot.background = element_blank(), 
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.grid.major.x = element_blank(),
                             panel.grid.major.y = element_blank())


hist<-ggplot(postmod2,aes(x=b_Prop.Change.Gn2))+ geom_density()+   transparent_theme
grob_hist<-ggplotGrob(hist)

bars + annotation_custom(grob = grob_hist,  xmin=min(postmod2$b_Prop.Change.Gn2), xmax=max(postmod2$b_Prop.Change.Gn2))



#how many negative, neutral, positive?? 
neutral<-subset(slopes, round(slopes$Estimate.Prop.Change.Gn2)==0)
 pos<-subset(slopes, round(slopes$Estimate.Prop.Change.Gn2)>0)
 neg<-subset(slopes, round(slopes$Estimate.Prop.Change.Gn2)<0)
 
 #can also look at CIs crossing zero or not (more accurate)
 pos<-subset(slopes, slopes$Q2.5.Prop.Change.Gn2>0 & slopes$Q97.5.Prop.Change.Gn2>0)
 #5643 or 25 % 
 neg<-subset(slopes, slopes$Q2.5.Prop.Change.Gn2<0 & slopes$Q97.5.Prop.Change.Gn2<0)
 #148 or <1% %
 sig<-rbind(pos, neg)
 neutral<-anti_join(slopes, sig) #17127 or 75%
 
 
 
 #random forest output----
 #ordinal forest
 conf_mat<-as.data.frame(conf_mat)
 
 ggplot(data=conf_mat, aes(y = predicted, x=true, fill=Freq)) + 
   geom_tile()+ 
   ylab("Predicted") + xlab("Observed") +
   theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

 ggplot(data=subset(conf_mat,true!=0& predicted!=0), aes(y = predicted, x=true, fill=Freq)) + 
   geom_tile()+ 
   ylab("Predicted") + xlab("Observed") +
   theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")
 