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
#try using suite of climate params from BGC projections (Ceres)
climrVars <- c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
               "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI")

#derived vars 
climPredictors <- c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
                    "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS")

var_names<-c(c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "Species") , climPredictors)
#"Tmax_sm", "TD", "PPT_sm", "DD5_sp") 

#create dataset with all preds for rf model 
tree_dat_sub<-select(tree_dat, var_names)
tree_dat_sub<-na.omit(tree_dat_sub)  #remove NAs

#which have nas?
nas<-anti_join(tree_dat, tree_dat_sub)
nas<-select(nas, var_names)

#check for correlation among climate vars- maybe not needed?
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

# Split the dataset into training and testing sets
set.seed(123) 
train_indices <- sample(1:nrow(tree_dat_sub), 0.7 * nrow(tree_dat_sub)) # 70% of the data for training
train_data <- tree_dat_sub[train_indices, ]
test_data <- tree_dat_sub[-train_indices, ]

#Fit the ordinal random forest model
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

#for non- zeroes- mis predicts 2s as 1s more often than correctly as 2s, predicts 3s as 0 or 2 as much as correctly as 3
ggplot(data=subset(conf_mat_df5,  predicted!=0 & true!=0), aes(y = predicted, x=true, fill=Freq)) + 
  geom_tile()+ 
  ylab("Predicted") + xlab("Observed") +
  theme_bw()+  scale_fill_gradient(low="white", high="darkgreen")

#look at feature importance 
importancedf<-as.data.frame(RFord$varimp)

#train model on ALL data -Final model
RFord <- ordfor(depvar = "cover_rank", mtry=14, data = tree_dat_sub, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodel8.Rdata")

#predict on all data (i.e. fitted values)
preds <- predict(RFord, newdata=tree_dat_sub)
fitted<-preds$ypred
preds$classprobs#what is this?

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

tree_dat_sub$fitted<-fitted




