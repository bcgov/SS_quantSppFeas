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
#load(file="data/tree_data_cleaned.Rdata") #update to OS path #13
load(file="data/tree_data_cleaned_wzeros.Rdata") #with absences #update to OS path #13
tree_dat<-tree_dat_wzeros #rename
rm(tree_dat_wzeros)
gc()

#climate variables----
#get all climate variables 
vars<-climr::variables #look up table for vars 
var_names<-vars$Code
climrVars <- c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
               "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI")
var_names<-c(c("TotalA", "Species", "NutrientRegime_clean", "MoistureRegime_clean", "Species") , climrVars)

#create dataset with all preds for rf model 
tree_dat_sub<-dplyr::select(tree_dat, var_names, Latitude, Longitude)

#tree_dat_sub<-na.omit(tree_dat_sub)  #remove NAs
gc()
#set Nutrient and Moisture Regimes as ordinal factor
str(tree_dat_sub)
tree_dat_sub$NutrientRegime_clean<-ordered(tree_dat_sub$NutrientRegime_clean, levels = c("A", "B", "C", "D", "E", "F"))
tree_dat_sub$MoistureRegime_clean<-ordered(tree_dat_sub$MoistureRegime_clean, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8"))
str(tree_dat_sub$MoistureRegime_clean)#good 
str(tree_dat_sub$NutrientRegime_clean)#good 

#response variables----
# Convert totalA to abundance classes
tree_dat_sub<-mutate(tree_dat_sub, TotalA_class=case_when(
  TotalA>50 & TotalA<=100~5,
  TotalA>25 & TotalA<=50~4,
  TotalA>=11 & TotalA<=25~3,
  TotalA>=5 & TotalA<=11~2,
  TotalA>0 & TotalA<5~1,
  #TotalA>0 & TotalA<1~1,
  TotalA==0 ~0, 
  TRUE~NA))


#look at distribution of classes 
hist(tree_dat_sub$TotalA_class) 
knitr::kable(group_by(tree_dat_sub, TotalA_class)%>%summarise(counts=n()))

#single edatope models----
#create a unique variable for the 15 edatopic spaces to train on 
tree_dat_sub$edatopex<-paste(tree_dat_sub$NutrientRegime_clean, tree_dat_sub$MoistureRegime_clean, sep="")
tree_dat_sub<- mutate(tree_dat_sub, edatope=case_when(edatopex=="C3"|edatopex=="C4"~"C34",
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
                     
                     
#check<-dplyr::select(tree_dat_sub, MoistureRegime_clean, NutrientRegime_clean, edatope, edatopex) 
unique(tree_dat_sub$edatope)


#save final modeling dataset
save(tree_dat_sub, file='data/final_model_data.Rdata')
rm(tree_dat)
gc()


#fit models 
library(ranger)

#remove continuous abundance for modeling 
tree_dat_sub$TotalA<-NULL 
#remove lat/long values for modeling 
tree_dat_sub$Latitude
tree_dat_sub$Longitude

# List of unique species & sites (edatopic spaces)
site_list<-unique(tree_dat_sub$edatope)
species_list <- unique(tree_dat_sub$Species)
output_dir<- "outputs/ranger/RFregression_classes/Single_edatope"

#just start with 1 spp & 3 sites to test that the loop works
#species_list<-species_list[14] 
#site_list<-site_list[c(1:2)] 

# Loop through each species, edatope and fit a model
for (site in site_list) {                                        
  for (species in species_list) {
  # Subset the data for the current species x edatope
  model_data <- tree_dat_sub %>% 
    filter(site == site, species == species) %>%
    dplyr::select(-edatope, -edatopex)#remove edatope as predictor variable (keep SMR, SNR) 
  
 # Fit RF model 
 RF<- ranger::ranger(TotalA_class~ .,data = model_data, mtry=20, classification = F) 
  
  # Create a filename for saving the model
  key<- paste(species, site, sep = "_")
  filename <- file.path(output_dir, paste0("RF_", key, ".RData"))
  
  # Save the model output to a file
  save(RF, file = filename)
  
 }
}


# Cross-validation
library(caret)
folds <- 5  # number of cross-validation folds

# Load data if not already loaded
#load(file='data/final_model_data.Rdata')
species_list <- unique(tree_dat_sub$Species)
site_list <- unique(tree_dat_sub$edatope)

# Remove parameters not in the model
tree_dat_sub <- select(tree_dat_sub, -Latitude, -Longitude)

# Initialize nested containers for models and predictions by species
container_model <- vector("list", length(species_list))
container_pred <- vector("list", length(species_list))

# Create a container for performance metrics
performance_metrics <- data.frame(Species = character(),
                                  Site = character(),
                                  RMSE = numeric(),
                                  R2 = numeric(),
                                  stringsAsFactors = FALSE)

# Loop through each species
for (species_index in seq_along(species_list)) {
  species <- species_list[species_index]
  
  # Filter data for the current species
  species_data <- tree_dat_sub %>% filter(Species == species)
  
  # Loop through each site
  for (site_index in seq_along(site_list)) {
    site <- site_list[site_index]
    
    # Filter data for the current site
    site_data <- species_data %>% filter(edatope == site)
    

    # Create folds 
    set.seed(123)  # for reproducibility
    if (nrow(site_data) < 5) next 
    cvIndex <- createFolds(1:nrow(site_data), k = 5, list = TRUE, returnTrain = TRUE)
    
    # Remove edatopes
    site_data <- select(site_data, -edatopex, -edatope)
    
    # Initialize containers for this species and site
    container_model[[species_index]][[site_index]] <- vector("list", length(cvIndex))
    container_pred[[species_index]][[site_index]] <- vector("list", length(cvIndex))
    
    # Iterate through the cross-validation folds
    for (i in 1:length(cvIndex)) {
      # Define training and evaluation data
      train_data <- site_data[cvIndex[[i]], ]
      eval_data <- site_data[-cvIndex[[i]], ]
      
            # Train the model
      rf <- ranger::ranger(TotalA_class ~ ., data = train_data, mtry = 20, classification = FALSE) 
      
      # Predict on the hold-out test set
      pred <- predict(rf, eval_data)$predictions
      
      # Store results in the appropriate species and site container
      container_model[[species_index]][[site_index]][[i]] <- rf
      container_pred[[species_index]][[site_index]][[i]] <- pred
      
      # Calculate true labels
      true <- eval_data$TotalA_class 
      
      # Calculate performance metrics
      rmse <- sqrt(mean((pred - true) ^ 2))  # Root Mean Squared Error
      r2 <- cor(pred, true) ^ 2  # R-squared
      
      # Store performance metrics
      performance_metrics <- rbind(performance_metrics, 
                                   data.frame(Species = species,
                                              Site = site,
                                              RMSE = rmse,
                                              R2 = r2,
                                              stringsAsFactors = FALSE))
    }
  }
}
performance_metrics_single<-group_by(performance_metrics,Species, Site)%>%summarise(RMSE=mean(RMSE), R2=mean(R2))


#multi edatope models----
#Now run model for each species with all edatopes 
species_list <- unique(tree_dat_sub$Species)
output_dir<- "outputs/ranger/RFregression_classes"

#just start with a few spp  to test that the loop works
#species_list<-species_list[12:16] 

#tree_dat_sub<-dplyr::select(tree_dat_sub, -edatope, -edatopex)#remove edatope as predictor variable (keep SMR, SNR) 

# Loop through each species, site and fit a model
  for (species in species_list) {
    # Subset the data for the current species 
    model_data <- tree_dat_sub %>% 
      filter(species == species) 
    
    # Fit RF model 
    RF<- ranger::ranger(TotalA_class~ .,data = model_data, mtry=20, classification = F) 
    
    # Create a filename for saving the model
    filename <- file.path(output_dir, paste0("RF_", species, ".RData"))
    
    # Save the model output to a file
    save(RF, file = filename)
    
  }

#Cross Validation
library(caret)
#cross validation, stratified on edatope to ensure that each group 
# is equally distributed over the cross-validation folds
folds <- 5 # for <nfold> cross-validation

#load data if not already
#load(file='data/final_model_data.Rdata')

species_list <- unique(tree_dat_sub$Species)

# Initialize nested containers for models and predictions by species
container_model <- vector("list", length(species_list))
container_pred <- vector("list", length(species_list))

# Create a container for performance metrics
performance_metrics <- data.frame(Species = character(),
                                  cor = numeric(),
                                  R2 = numeric(),
                                  stringsAsFactors = FALSE)

#remove params not in model 
tree_dat_sub<-select(tree_dat_sub, -Latitude, -Longitude, -TotalA)

# Loop through each species
for (species_index in seq_along(species_list)) {
  species <- species_list[species_index]
  
  # Filter data for the current species
  species_data <- tree_dat_sub %>% filter(Species == !!species)
  
  # Create stratified folds based on 'edatope' for the current species
  set.seed(123)  # for reproducibility
  cvIndex <- createFolds(species_data$edatope, k = folds, list = TRUE, returnTrain = TRUE)
  
  #remove edatope from dataset
  species_data<-select(species_data, -edatopex, -edatope)
  
  # Initialize containers for this species
  container_model[[species_index]] <- vector("list", length(cvIndex))
  container_pred[[species_index]] <- vector("list", length(cvIndex))
  
  # Iterate through the cross-validation folds
  for (i in 1:length(cvIndex)) {
    
    # Define training and evaluation data
    train_data <- species_data[cvIndex[[i]], ]
    eval_data <- species_data[-cvIndex[[i]], ]
    
    # Ensure there is enough data in the evaluation set
    if (nrow(eval_data) == 0) next
    
    # Train the model
    rf <- ranger::ranger(TotalA_class~ .,data = train_data, mtry=20, classification = F) 
    
    # Predict on the hold-out test set
    pred <- predict(rf, eval_data)$predictions
    
    # Store results in the appropriate species container
    container_model[[species_index]][[i]] <- rf
    container_pred[[species_index]][[i]] <- pred
    
    # Calculate true labels
    true <- eval_data$TotalA_class 
    
    # Calculate performance metrics
    rmse <- sqrt(mean((pred - true) ^ 2))  # Root Mean Squared Error
    r2 <- cor(pred, true) ^ 2  # R-squared
    
    # Store performance metrics
    performance_metrics <- rbind(performance_metrics, 
                                 data.frame(Species = species,
                                            RMSE = rmse,
                                            R2 = r2,
                                            stringsAsFactors = FALSE))
  }
}

performance_metrics<-group_by(performance_metrics,Species)%>%summarise(RMSE=mean(RMSE), R2=mean(R2))



#combine performance_metrics single and multi edatope models 
performance_metrics$Site<-"all"
performance_metrics<-relocate(performance_metrics, Site, .after = Species)
performance_metrics<-rbind(performance_metrics, performance_metrics_single)
