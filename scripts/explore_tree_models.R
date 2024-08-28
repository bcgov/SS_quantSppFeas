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
load(file="data/tree_data_cleaned_wzeros.Rdata") #with absences #update to local path 
tree_dat<-tree_dat_wzeros #rename
rm(tree_dat_wzeros)
gc()

#plot climate variables----

#get all climate variables 
vars<-climr::variables #look up table for vars 
var_names<-vars$Code
#add response variable and a few other preds to all climate data 
#var_names<-c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "SlopeGradient", "Aspect", "year", "bgc" , "Elevation" ,var_names)

#subset predictor vars
var_names<-c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean","Tmax_sm", "TD", "PPT_sm", "DD5_sp")

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
RFord <- ordfor(depvar = "cover_rank", mtry=7, data = train_data, ntreefinal = 1000)  
save(RFord, file="outputs/ordinalForest/RFordmodel5.Rdata")

# Predict on the test set
preds <- predict(RFord, newdata=test_data)

#confusion matrix 
conf_mat<-table('true'=test_data$cover_rank, 'predicted'=preds$ypred) 
accuracy <- sum(diag(conf_mat)) / sum(conf_mat) #acc = 0.88

