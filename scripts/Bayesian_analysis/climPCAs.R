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

library(vegan)
library(ape)

#load feas tables with climate data 
load(file="data/feas_abund_clim_data.Rdata")

#remove NAs for ordination
feas.dat.clim<-na.omit(feas.dat.clim)

#removing large climate outlier plot in Mhmm1/101a
feas.dat.clim<-subset(feas.dat.clim, PlotNumber!= "S95CT71")


#select climate variables from BGC model
climrVars = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
              "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
              "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
              "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", "PPT_JAS", "CMD.total")

climdat<- select(feas.dat.clim, climrVars) %>%
  #mutate(PASlog=log(PAS), PAS_splog=log(PAS_sp),PPT_MJlog=log(PPT_MJ), PPT_JASlog=log(PPT_JAS))%>% #log precip vars
  #select(-PAS, -PAS_sp, -PPT_MJ, -PPT_JAS) %>%
  scale(.) # standardize the data

# Perform principal component analysis (PCA) on the standardized data
climpca<-  prcomp(climdat)

# Print a summary of the PCA results
#summary(climpca)
#biplot(climpca)
#scores(climpca)

# Convert the PCA results to a tibble and add a column for observation names
climpca_df <- as_tibble(climpca$x) %>% 
  mutate(observation = rownames(.))%>%select(PC1, PC2, PC3, observation)

#merge PCs back with feas data---- 
feas.dat.clim<-add_rownames(feas.dat.clim, "observation")
feas.dat.clim<-left_join(climpca_df, feas.dat.clim)

#save output
save(feas.dat.clim, file="data/feas_abund_clim_data.Rdata")


#PCs in in vegan
#climpca2 <- vegan::rda(climdat, scale=F) #already scaled above 
#scores <- scores(climpca2)
#biplot(climpca) 
#eigenvals(climpca2)

#loading scores of climate params onto PCA axes 
#climloads<-envfit(climpca, climdat, permutations = 999, strata = NULL)
#climloads
#scrs <- as.data.frame(scores(climloads, display = "vectors"))
#scrs<-add_rownames(scrs, "clim")
#write.csv(scrs, "ClimPCA_loadingscores.csv")

#ggplot(scrs) +
#  geom_point(mapping = aes(x = PC1, y = PC2)) +
#  coord_fixed() + ## need aspect ratio of 1!
#  geom_segment(data = scrs,
#               aes(x = 0, xend = PC1, y = 0, yend = PC2),
#               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#  geom_text(data = scrs, aes(x = PC1, y = PC2, label = clim),
#            size = 3)


