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
tree_dat_sub<-dplyr::select(tree_dat, var_names)
#tree_dat_sub<-na.omit(tree_dat_sub)  #remove NAs
gc()

#response variables----
# Convert totalA to abundance classes
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


#train random forest models by spp----
#remove continuous abundance for modeling 
tree_dat_sub$TotalA<-NULL 

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
                                                #edatopex=="F5"|edatopex=="F6"~"F56", #F not in ByBEC
                                                #                  edatopex=="F3"~"F3", #F not in BYBEC
                                                                            TRUE~NA))
                     
                     
check<-dplyr::select(tree_dat_sub, MoistureRegime_clean, NutrientRegime_clean, edatope, edatopex) 
unique(tree_dat_sub$edatope)

#fit models 
library(ranger)

# List of unique species & sites (edatopic spaces)
site_list<-unique(tree_dat_sub$edatope)
species_list <- unique(tree_dat_sub$Species)
output_dir<- "outputs/ranger/RFregression_classes"

#just start with 1 spp & 3 sites to test that the loop works
#species_list<-species_list[14:16] 
#site_list<-site_list[c(1:2,5)] 

# Loop through each species, site and fit a model
for (site in site_list) {                                        
  for (species in species_list) {
  # Subset the data for the current species x edatope
  model_data <- tree_dat_sub %>% 
    filter(site == site, species == species) %>%
    dplyr::select(-edatope, -edatopex)#remove edatope as predictor variable (keep SMR, SNR) 
  
 # Fit RF model 
 RF<- ranger::ranger(TotalA_class~ .,data = model_data, mtry=20, classification = F) #update mtry based on which climrVars
  
  # Create a filename for saving the model
  key<- paste(species, site, sep = "_")
  filename <- file.path(output_dir, paste0("RF_", key, ".RData"))
  
  # Save the model output to a file
  save(RF, file = filename)
  
 }
}


#Now run model for each species with all edatopes 
species_list <- unique(tree_dat_sub$Species)
output_dir<- "outputs/ranger/RFregression_classes"

#just start with 3 spp  to test that the loop works
species_list<-species_list[15:17] 

tree_dat_sub<-dplyr::select(tree_dat_sub, -edatope, -edatopex)#remove edatope as predictor variable (keep SMR, SNR) 

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

