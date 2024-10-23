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
#load(file="data/tree_data_cleaned.Rdata") #cover data only 
load(file="data/tree_data_cleaned_wzeros.Rdata") #with absences #update to OS path #13
tree_dat<-tree_dat_wzeros #rename
rm(tree_dat_wzeros)
gc()

#climate variables----
#get all climate variables 
vars<-climr::variables #look up table for vars 
var_names<-vars$Code
#add response variable and a few other preds to all climate data 
#var_names<-c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "SlopeGradient", "Aspect", "year", "bgc" , "Elevation" ,var_names)

#subset predictor vars
#try using suite of 17 climate params from BGC projections (Ceres)
climrVars <- c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
               "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI")

#derived vars (11)
climPredictors <- c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
                    "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS")

#tmin, tmax, and precip for the four seasons- for diagnosis
climSimple<-c("PPT_at", "PPT_wt", "PPT_sp", "PPT_sm", 'Tmax_at', "Tmax_wt", "Tmax_sp", "Tmax_sm",
              'Tmin_at', "Tmin_wt", "Tmin_sp", "Tmin_sm")
              
#choose one 
var_names<-c(c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "Species") , climrVars)
#var_names<-c(c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "Species") , climPredictors)
#var_names<-c(c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "Species") , climSimple)

#create dataset with all preds for rf model 
tree_dat_sub<-dplyr::select(tree_dat, var_names)

#check for correlation among climate vars----
#maybe not needed?
#kiri suggests throw out r>0.9
#Colin- correlation at entire dataset scale not relevant due to splitting by RF algorithm
#pairs(tree_dat_sub[,c(5:15)])

cor.test(tree_dat_sub$PPT_05, tree_dat_sub$PPT_06)#correlated ~0.7 remove
cor.test(tree_dat$DD5_sp, tree_dat$Tmax_sm)#correlated ~0.65 better
plot(tree_dat$CMI, tree_dat$TD)
cor.test(tree_dat$NFFD, tree_dat$TD) #r~ 0.45 ok 

plot(tree_dat$DD5, tree_dat$Tmax_sm)
cor.test(tree_dat$DD5, tree_dat$Tmax_sm)#correlated ~0.7 remove
cor.test(tree_dat$DD5_sp, tree_dat$Tmax_sm)#correlated ~0.65 better
plot(tree_dat$CMI, tree_dat$TD)
cor.test(tree_dat$NFFD, tree_dat$TD) #r~ 0.45 ok 

cor.test(tree_dat_sub$PPT_at, tree_dat_sub$PPT_wt) #0.97 too high
cor.test(tree_dat_sub$CMI, tree_dat_sub$PPT_wt) #0.98 
cor.test(tree_dat_sub$CMI, tree_dat_sub$PPT_at) #0.99
#remove ppt wt & ppt at 

cor.test(tree_dat_sub$CMD, tree_dat_sub$CMD_07) #0.94
cor.test(tree_dat_sub$PPT_05, tree_dat_sub$PPT_07) #0.8 ok
cor.test(tree_dat_sub$PPT_06, tree_dat_sub$PPT_07) #0.88 ok
cor.test(tree_dat_sub$PPT_06, tree_dat_sub$PPT_05) #0.9 
cor.test(tree_dat_sub$PPT_08, tree_dat_sub$PPT_07) #0.9 
cor.test(tree_dat_sub$PPT_08, tree_dat_sub$PPT_09) #0.93 
cor.test(tree_dat_sub$PPT_05, tree_dat_sub$PPT_09) #0.93 

cor.test(tree_dat_sub$DDsub0_at, tree_dat_sub$DDsub0_wt) #0.96 

#11 derived vars-test correlations 
cor.test(tree_dat_sub$PPT_MJ, tree_dat_sub$PPT_JAS) #0.92 too high?
cor.test(tree_dat_sub$CMI, tree_dat_sub$PPT_JAS) #0.92
cor.test(tree_dat_sub$CMD.total, tree_dat_sub$AHM) #0.96


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
  mtry=9, 
  scale.permutation.importance = TRUE)

save(RFmodel_allspp, file="outputs/ranger/RFmodel7.Rdata")
print(RFmodel_allspp) #R2 max=0.3

# Predict on the test set
preds <- predict(RFmodel_allspp, data = test_data)

#regress
plot(test_data$TotalA, preds$predictions)
cor.test(test_data$TotalA, preds$predictions)#r=0.52, r2=0.27 (so about the same as OOB)


#RF regression ----
#start with one species & edatope and 6 classes   
tree_dat_sub<-mutate(tree_dat_sub, TotalA_class=case_when(
                                         TotalA>50 & TotalA<=100~5,
                                         TotalA>25 & TotalA<=50~4,
                                         TotalA>=5 & TotalA<=25~3,
                                         TotalA>=1 & TotalA<5~2,
                                         TotalA>0 & TotalA<1~1,
                                           TotalA==0 ~0, 
                                           TRUE~NA))

#look at distribution of classes 
hist(tree_dat_sub$TotalA_class) 
knitr::kable(group_by(tree_dat_sub, TotalA_class)%>%summarise(counts=n()))


#run for one spp/edatope space 
tree_dat_sub_Hw<-subset(tree_dat_sub, Species=="TSUGHET")
#now try predicting just one edatope
tree_dat_sub_Hwzonal<-subset(tree_dat_sub_Hw, NutrientRegime_clean=="C")%>%subset(MoistureRegime_clean=="4"|MoistureRegime_clean=="3")
#remove continuous abundance for modeling 
tree_dat_sub_Hwzonal$TotalA<-NULL 

# Split the dataset into training and testing sets
set.seed(123) 
train_indices <- sample(1:nrow(tree_dat_sub_Hwzonal), 0.7 * nrow(tree_dat_sub_Hwzonal)) # 70% of the data for training
train_data <- tree_dat_sub_Hwzonal[train_indices, ]
test_data <- tree_dat_sub_Hwzonal[-train_indices, ]

#run the model - 12 simple climate params 
RFmodel_HWzonal <- ranger::ranger(
  TotalA_class~ .,
  data = train_data,
  #splitrule = "maxstat",  #recommended for accurate importance feature ranking
  #importance = 'permutation',
  mtry=15, classification = F) 
  #scale.permutation.importance = TRUE)

# Predict on the test set
preds <- predict(RFmodel_HWzonal, data = test_data)

#regress
plot(test_data$TotalA_class, preds$predictions)
cor.test(test_data$TotalA_class, preds$predictions)#r=0.76 %

#try the same model but with classification 
RFmodel_HWzonal_class <- ranger::ranger(
  TotalA_class~ .,
  data = train_data,
  #splitrule = "maxstat",  #recommended for accurate importance feature ranking
  #importance = 'permutation',
  mtry=15, classification = T) 
#scale.permutation.importance = TRUE)
# Predict on the test set
preds2 <- predict(RFmodel_HWzonal_class, data = test_data)

#confusion matrix 
conf_mat<-table('true'=test_data$TotalA_class, 'predicted'=preds2$predictions) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.66 
conf_mat

#Does better with regression


#run the model - 11 derived climate params 
RFmodel_HWzonal2 <- ranger::ranger(
  TotalA_class~ .,
  data = train_data,
  #splitrule = "maxstat",  #recommended for accurate importance feature ranking
  #importance = 'permutation',
  mtry=14, classification = F) 

# Predict on the test set
preds <- predict(RFmodel_HWzonal2, data = test_data)

#regress
plot(test_data$TotalA_class, preds$predictions)
cor.test(test_data$TotalA_class, preds$predictions)#r=0.815

#run the model - 17 climate params 
RFmodel_HWzonal3 <- ranger::ranger(
  TotalA_class~ .,
  data = train_data,
  #splitrule = "maxstat",  #recommended for accurate importance feature ranking
  #importance = 'permutation',
  mtry=20, classification = F) 

# Predict on the test set
preds <- predict(RFmodel_HWzonal3, data = test_data)

#regress
plot(test_data$TotalA_class, preds$predictions)
cor.test(test_data$TotalA_class, preds$predictions)#r=0.82

#compare modsa
print(RFmodel_HWzonal3) #R2=0.64
print(RFmodel_HWzonal2) #R2=0.63
print(RFmodel_HWzonal) #R2=0.64


#run the model on all data for training- 17 climate params 
RFmodel_HWzonal3 <- ranger::ranger(
  TotalA_class~ .,
  data = tree_dat_sub_Hwzonal,
  #splitrule = "maxstat",  #recommended for accurate importance feature ranking
  #importance = 'permutation',
  mtry=20, classification = F) 

save(RFmodel_HWzonal3, file="outputs/ranger/RFmodelHwzonal.Rdata")

print(RFmodel_HWzonal3)#R2=0.68

#ordinal forest models----
# Convert totalA to an ordered factor
#breaks from plot data (explore_veg_data.R)
#newfeas  meds  mean
#       1    17 21.8 
#       2    10 16.2 
#       3     5  9.64
#       4     2  2   
#       5     2  4.67
rm(tree_dat)
gc()
tree_dat_sub<-mutate(tree_dat_sub, cover_rank=case_when(TotalA>20~1, 
                                                        TotalA>0 & TotalA<7.5~3,
                                                        TotalA==0 ~0, 
                                                        TRUE~2))
hist(tree_dat_sub$cover_rank)

group_by(tree_dat_sub, cover_rank)%>%summarise(counts=n())#about even for 1, 2,3

#set as ordinal factor
tree_dat_sub$cover_rank<-ordered(tree_dat_sub$cover_rank, levels = c(0, 3, 2, 1))
str(tree_dat_sub$cover_rank)#good 

#remove continuous response
tree_dat_sub$TotalA<-NULL


save(tree_dat_sub,file="data/ordforest_data.Rdata")#17 climate vars
save(tree_dat_sub,file="data/ordforest_data2.Rdata")#11 derived climate vars 

# Split the dataset into training and testing sets
set.seed(123) 
train_indices <- sample(1:nrow(tree_dat_sub), 0.7 * nrow(tree_dat_sub)) # 70% of the data for training
train_data <- tree_dat_sub[train_indices, ]
test_data <- tree_dat_sub[-train_indices, ]

#Fit the ordinal random forest model----
#load datasets
load(file="data/ordforest_data.Rdata")#17 climate vars
load(file="data/ordforest_data2.Rdata")#11 derived climate vars 

library(ordinalForest)
RFord <- ordfor(depvar = "cover_rank", mtry=14, data = train_data, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodel7.Rdata")

# Predict on the test set
preds5 <- predict(RFord5, newdata=test_data)

#confusion matrix 
conf_mat<-table('true'=test_data$cover_rank, 'predicted'=preds$ypred) 
accuracy5 <- sum(diag(conf_mat5)) / sum(conf_mat5) #acc = 0.88
conf_mat

#plot CM for all data
conf_mat_df<-as.data.frame(conf_mat)
ggplot(data=conf_mat_df5, aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

ggplot(data=subset(conf_mat_df5,  predicted!=0 & true!=0), aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#look at feature importance 
importancedf<-as.data.frame(RFord$varimp)

#train model on ALL data 
RFord <- ordfor(depvar = "cover_rank", mtry=20, data = tree_dat_sub, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodel9.Rdata")

#predict on all data (i.e. fitted values)
preds <- predict(RFord, newdata=tree_dat_sub)
fitted<-preds$ypred
#preds$classprobs#what is this?

#confusion matrix 
conf_mat<-table('true'=tree_dat_sub$cover_rank, 'predicted'=preds$ypred) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.99
conf_mat #looks good 328 total incorrect (~0.2%)- mostly 2s being predicted as 1s or 3s

#plot CM for all data
conf_mat_df<-as.data.frame(conf_mat)
ggplot(data=conf_mat_df, aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#for non- zeroes- mis predicts 2s as 1s
conf_mat_df_sub<-subset(conf_mat_df, predicted!=0 & true!=0)

ggplot(data=conf_mat_df_sub, aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")


#try fitting ordinal forest for one species & zonal (Hw)
tree_dat_sub_Hwzonal<-subset(tree_dat_sub, Species=="TSUGHET"& NutrientRegime_clean=="C")%>%subset(MoistureRegime_clean=="4"|MoistureRegime_clean=="3")
#tree_dat_sub_Hwzonal<-dplyr::select(tree_dat_sub_Hwzonal, -Species, -NutrientRegime_clean, -MoistureRegime_clean) #remove categorical vars

RFord <- ordfor(depvar = "cover_rank", mtry=14, data = tree_dat_sub_Hwzonal, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodelHwZonal.Rdata")

#predict on i.e. fitted values
preds <- predict(RFord, newdata=tree_dat_sub_Hwzonal)
fitted<-preds$ypred
#confusion matrix 
conf_mat<-table('true'=tree_dat_sub_Hwzonal$cover_rank, 'predicted'=preds$ypred) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.997
conf_mat 

#plot CM for all data
conf_mat_df<-as.data.frame(conf_mat)
ggplot(data=conf_mat_df, aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#for non- zeroes- mis predicts 2s as 1s
conf_mat_df_sub<-subset(conf_mat_df, predicted!=0 & true!=0)

ggplot(data=conf_mat_df_sub, aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#try fitting ordinal forest for one spp (all edatopes)
tree_dat_sub_Hw<-subset(tree_dat_sub, Species=="TSUGHET")
tree_dat_sub_Hw$Species<-NULL
RFordHw2 <- ordfor(depvar = "cover_rank", mtry=13, data = tree_dat_sub_Hw, ntreefinal = 1000)  
save(RFordHw2, file="outputs/ordinalForest/RFordmodelHw2.Rdata")

preds <- predict(RFordHw2, newdata=tree_dat_sub_Hw)
predszonal <- predict(RFordHw2, newdata=tree_dat_sub_Hwzonal)

#confusion matrix 
conf_mat<-table('true'=tree_dat_sub_Hw$cover_rank, 'predicted'=preds$ypred) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.997

#now try predicting just one edatope
tree_dat_sub_Hwzonal<-subset(tree_dat_sub_Hw, NutrientRegime_clean=="C")%>%subset(MoistureRegime_clean=="4"|MoistureRegime_clean=="3")
preds <- predict(RFordHw2, newdata=tree_dat_sub_Hwzonal)
#confusion matrix 
conf_mat<-table('true'=tree_dat_sub_Hwzonal$cover_rank, 'predicted'=predszonal$ypred) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.71
conf_mat



