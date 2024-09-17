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
library(terra)
library(data.table)
library(climr)
library(sf)
library(raster)#masks dplyr select!!

#Use PRISM DEM to call in downscaled climr data----
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="") #must have VPN connected- 800 m DEM
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))

## convert the DEM to a data.frame
my_grid <- as.data.frame(dem.bc, cells = TRUE, xy = TRUE)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

## climr call- This will return the observed 1961-1990 climates for the raster grid points. 
clim.bcv <- downscale(
  xyz = my_grid, 
  out_spatial = TRUE,
  #obs_periods = "2001_2020", 
  vars = c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
               "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", 
                "CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp"))

saveRDS(clim.bcv, file="data/spatial/clim.bcv.RData")
#clim.bcv<-readRDS(file="data/spatial/clim.bcv.RData")
plot(clim.bcv)
gc()

#rasterize from vector 
my_rast <- rast(dem.bc) # use the DEM as a template raster
vars = c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
         "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", 
         "CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim.bcr<-rasterize(clim.bcv, my_rast, field= vars) #set background=0 for all cells that are NAs 
saveRDS(clim.bcr, file="data/spatial/clim.bcr.RData")
#clim.bcr<-readRDS(file="data/spatial/clim.bcr.RData")
plot(clim.bcr$DD5)#looks good 
gc()

#aggregate up to 2km resolution
clim.bcr2k <- aggregate(clim.bcr, fact=3, fun=mean)
saveRDS(clim.bcr2k, file="data/spatial/clim.bcr2k.RData")
#clim.bcr2k<-readRDS(file="data/spatial/clim.bcr2k.RData")
plot(clim.bcr2k$DD5)#looks good - can't really see resolution change
gc()

#clip to BC boundary
bcboundary<-bcmaps::bc_bound_hres()
#st_crs(bcboundary)
plot(st_geometry(bcboundary))
bcbound.reproj <- st_transform(bcboundary, st_crs(4326)) #reproject to wgs84  
clim.bcr2kcrop<-terra::crop(x = clim.bcr2k, y = bcbound.reproj) #crop climr output to BC boundary 

plot(clim.bcr2kcrop$DD5)#looks good- cropped slightly
gc()


#pull back out to df
clim.bc <- as.data.frame(clim.bcr2kcrop, cells = TRUE, xy = TRUE, na.rm=F)
clim.bc<-rename(clim.bc, id=cell, lon=x, lat=y)

#where are the NAs?
#clim.bc <- as.data.frame(clim.bcr2kcrop, cells = TRUE, xy = TRUE)
#clim.bcwna <- as.data.frame(clim.bcr2kcrop, cells = TRUE, xy = TRUE, na.rm=F)
#climnas<-anti_join(clim.bcwna, clim.bc) #129284 NAs 
#climnas<-rename(climnas, id=cell, lon=x, lat=y)%>%dplyr::select(id, lat, lon)
#plot(clim.bcr2kcrop$DD5)
#points(x = climnas$lon, 
#        y = climnas$lat, 
#        col = "red", 
#        cex = 0.5)
#write.csv(climnas, 'data/gpsnas.csv')

#additional vars needed -"DD_delayed", "PPT_MJ", "PPT_JAS", "CMD.total", "CMDMax"----
clim.bc<-as.data.table(clim.bc) #requires data table 
source("scripts/addVars.R") 
addVars(clim.bc) #why not showing additional variables in environment??
save(clim.bc, file="data/clim.bc2k.RData")

#predict over full climate surface with trained OF model(s)----
#load climate all BC (from PRISM DEM) 
#load(file="data/clim.bc.RData") #800m 
load(file="data/clim.bc2k.RData") #2km
#fill all NAs with zero b/c predict function requires 
clim.bc[is.na(clim.bc)] <- 0

#load trained models 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")
RFord5<-RFord
load(file="outputs/ordinalForest/RFordmodel6.Rdata")
RFord6<-RFord
load(file="outputs/ordinalForest/RFordmodel7.Rdata")
RFord7<-RFord
rm(RFord)

#add edaphic info - start with Zonal for all 
clim.bc$NutrientRegime_clean<-"C" #zonal  
clim.bc$MoistureRegime_clean<-"4" #zonal  

#add species - start with Hw
clim.bc$Species<-"TSUGHET"

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim5<-dplyr::select(clim.bc, clim5)  
clim6<-c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
"PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", "CMI", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim6<-dplyr::select(clim.bc, clim6)  
clim7<-c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
        "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS", "NutrientRegime_clean", "MoistureRegime_clean", "Species")
clim7<-dplyr::select(clim.bc, clim7)  
gc()

#make predictions for all BC 
library(ordinalForest)
predsm5c4Hw <-predict(object = RFord5, newdata = clim5)
save(predsm5c4Hw, file="outputs/ordinalForest/preds/predsm5c4Hw.RData")
predsm6c4Hw <-predict(object = RFord6, newdata = clim6)
save(predsm6c4Hw, file="outputs/ordinalForest/preds/predsm6c4Hw.RData")
predsm7c4Hw <-predict(object = RFord7, newdata = clim7)
save(predsm7c4Hw, file="outputs/ordinalForest/preds/predsm7c4Hw.RData")

#make into df
predsm5c4Hw<-as.data.frame(predsm5c4Hw$ypred) 
predsm5c4Hw$id<-as.numeric(as.character(row.names(predsm5c4Hw)))
colnames(predsm5c4Hw)<- c("ypred", "id")

predsm6c4Hw<-as.data.frame(predsm6c4Hw$ypred) 
predsm6c4Hw$id<-as.numeric(as.character(row.names(predsm6c4Hw)))
colnames(predsm6c4Hw)<- c("ypred", "id")

predsm7c4Hw<-as.data.frame(predsm7c4Hw$ypred) 
predsm7c4Hw$id<-as.numeric(as.character(row.names(predsm7c4Hw)))
colnames(predsm7c4Hw)<- c("ypred", "id")

##add lat, long back in
predsm5c4Hw<-left_join(predsm5c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm6c4Hw<-left_join(predsm6c4Hw, dplyr::select(clim.bc, id, lat, lon))
predsm7c4Hw<-left_join(predsm7c4Hw, dplyr::select(clim.bc, id, lat, lon))

climnas<-read.csv("data/gpsnas.csv")
predsm5c4Hw<-filter(predsm5c4Hw, !id %in% climnas$id) # remove cells in ocean 
predsm6c4Hw<-filter(predsm6c4Hw, !id %in% climnas$id) # remove cells in ocean 
predsm7c4Hw<-filter(predsm7c4Hw, !id %in% climnas$id) # remove cells in ocean 


#plot predictions----
#BEC plot locations
load(file="data/tree_data_cleaned.Rdata")
Hw<-subset(tree_dat, Species=='TSUGHET') 

library(tidyterra)
#NOT WORKING
ggplot() +
  geom_spatraster(data = clim.bcr2kcrop, aes(fill = Tmax_sm))+
  scale_fill_whitebox_c(
    palette = "muted",
    na.value = "white") + 
    geom_point(data=predsm5c4Hw, aes(x=lon, y=lat, fill=ypred))

#base R
plot(clim.bcr2kcrop$DD5, legend=F, main = "Mod5 Hw c4 2km predited feasibility") #cropped raster
#overlay pred values
points(x = predsm5c4Hw$lon, predsm5c4Hw$lat,
       col=c("lightgrey","blue", "yellow", "green")[predsm5c4Hw$ypred],
       cex = 0.8)
#plot BEC plots
points(x = Hw$Longitude, Hw$Latitude,
       col='grey', pch=4,
       cex = 0.3)
legend("topright",   
       legend = levels(predsm5c4Hw$ypred), 
       # the legend
       fill =  c("lightgrey","blue", "yellow", "green"), 
       bty="n",cex=0.8, 
       inset=c(0.025, 0.1))
       
#mod 6
#load trained model 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-dplyr::select(clim.bc, clim5)  
clim5<-na.omit(clim5)  #remove NAs

#add edaphic info 
clim5$NutrientRegime_clean<-"C" #zonal  
clim5$MoistureRegime_clean<-"4" #zonal  
clim5c4<-clim5
clim5$NutrientRegime_clean<-"B" #dry poor 
clim5$MoistureRegime_clean<-"2" #dry poor
clim5b2<-clim5
clim5$NutrientRegime_clean<-"E" #wet rich 
clim5$MoistureRegime_clean<-"6" #wet rich 
clim5e6<-clim5

#add species 
clim5c4$Species<-"TSUGHET"
clim5b2$Species<-"TSUGHET"
clim5e6$Species<-"TSUGHET"

#make predictions for all BC 
rm(clim.bc)
gc()
predsm5c4Hw <- predict(object = RFord,newdata = clim5c4)
predsm5b2Hw <- predict(object = RFord,newdata = clim5b2)
predsm5e6Hw <- predict(object = RFord,newdata = clim5e6)



## populate the raster grid with prediction values
preds_map <- rast(dem.bc) # use the DEM as a template raster
preds_map[predsm5c4$id]] <- ds_out[PERIOD == "2001_2020", MAT] # populate the raster cells with the 2001-2020 mean annual temperature (MAT) values, using the `id` field as the link. 

plot(my_clim, main = "2001-2020 mean annual temperature (MAT)")


#mod 7
#load trained model 
load(file="outputs/ordinalForest/RFordmodel5.Rdata")

#subset to correct climate params 
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-dplyr::select(clim.bc, clim5)  
clim5<-na.omit(clim5)  #remove NAs

#add edaphic info 
clim5$NutrientRegime_clean<-"C" #zonal  
clim5$MoistureRegime_clean<-"4" #zonal  
clim5c4<-clim5
clim5$NutrientRegime_clean<-"B" #dry poor 
clim5$MoistureRegime_clean<-"2" #dry poor
clim5b2<-clim5
clim5$NutrientRegime_clean<-"E" #wet rich 
clim5$MoistureRegime_clean<-"6" #wet rich 
clim5e6<-clim5

#add species 
clim5c4$Species<-"TSUGHET"
clim5b2$Species<-"TSUGHET"
clim5e6$Species<-"TSUGHET"

#make predictions for all BC 
rm(clim.bc)
gc()
predsm5c4Hw <- predict(object = RFord,newdata = clim5c4)
predsm5b2Hw <- predict(object = RFord,newdata = clim5c4)
predsm5e6Hw <- predict(object = RFord,newdata = clim5c4)

## populate the raster grid with prediction values
preds_map <- rast(dem.bc) # use the DEM as a template raster
preds_map[predsm5c4$id]] <- ds_out[PERIOD == "2001_2020", MAT] # populate the raster cells with the 2001-2020 mean annual temperature (MAT) values, using the `id` field as the link. 

plot(my_clim, main = "2001-2020 mean annual temperature (MAT)")












clim.bc.stack<-raster::stack(clim.bcr) #stack 

#save(clim.bc.stack, file="data/clim.bc.stack.RData")
plot(clim.bcr) #looks good 

#add nutrient & moisture regimes 
a<-raster(ncol=3241, nrow=1680)
a<-setValues(a, "A")
clim.bc.stack<-addLayer(clim.bc.stack, a)
b<-rast(ncol=3241, nrow=1680)
c<-rast(ncol=3241, nrow=1680)
d<-rast(ncol=3241, nrow=1680)
e<-rast(ncol=3241, nrow=1680)
f<-rast(ncol=3241, nrow=1680)

#predict onto full climate surface using trained OF models 
#mod5
load(file="data/spatial/clim.bc.stack.RData")



load(file="outputs/ordinalForest/RFordmodel5.Rdata")
clim5<-c("Tmax_sm", "TD", "PPT_sm", "DD5_sp")
clim5<-select(clim.bc, clim5)  
clim5$NutrientRegime_clean<-"C" #use zonal here 
clim5$MoistureRegime_clean<-"3" #use zonal here 


spat_pred5<-raster::predict(clim5r, RFord, type="response")

rm(RFord)





---------------------
#add gps coords back in 
clim.bc<-left_join(clim.bc, my_grid)
clim.bc<-na.omit(clim.bc)  #remove NAs

clim.bc<-raster::rasterFromXYZ(clim.bc)

#lower resolution of latitute and longitude to 2 sig digits ~1km 
clim.bc$latlr<-round(clim.bc$lat, 3)
clim.bc$lonlr<-round(clim.bc$lon, 3)
#clim.bc<-group_by(clim.bc, latlr, lonlr)%>% mutate(elevlr=mean(elev)) #take mean elevs by GPS 

clim.bcx<-group_by(clim.bc, latlr, lonlr, elev)%>%   
  summarize(DD5_sp= mean(DD5_sp), Tmax_sm= mean(Tmax_sm), TD= mean(TD), PPT_sm=mean(PPT_sm))

#DDsub0_at= mean(DDsub0_at)
, "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", "PPT_08",
"PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", "SHM", "AHM", "NFFD", "PAS", 
"CMI", "Tmax_sm", "TD", "PPT_sm", "DD5_sp"))
gc()



#bring in feas tables---- 
#feas_tab<-read.csv("data/FeasibilityUpdates.csv")#downloaded from ByBEC 6/3/24
feas_tab<-read.csv("data/Feasibility_v12_15.csv")#most recent version from ccisr/data-raw/data_tables/


#update species naming in feas tables----
unique(sort(feas_tab$spp))

feas_tab<-mutate(feas_tab, Species= case_when(spp=="Ba"~"ABIEAMA",
                                                spp=="Ra" ~"ARBUMEN",
                                                spp=="Bg"~"ABIEGRA", 
                                                spp=="Bl"~"ABIELAS", 
                                                spp=="Mb"~"ACERMAC", 
                                                spp=="Dr"~"ALNURUB",
                                                spp=="Ra"~"ARBUMEN",  
                                                spp=="Ep"~"BETUPAP", 
                                                spp=="Yc"~"CALLNOO", 
                                                spp=="Lw"~"LARIOCC",  
                                                spp=="Se"~"PICEENE",
                                                spp=="Sw"~"PICEGLA",
                                                spp=="Sb"~"PICEMAR",
                                                spp=="Ss"~"PICESIT",
                                                spp=="Sxl"~"PICEXLU",
                                                spp=="Pa"~"PINUALB",
                                                spp=="Pl"~"PINUCON",
                                                spp=="Plc"~"PINUCON", #coastal
                                                spp=="Pli"~"PINUCON", #interior
                                                spp=="Pyc"~"PINUPON",
                                                spp=="Pyi"~"PINUPON",
                                                spp=="Pw"~"PINUMON",
                                                spp=="Acb"~"POPUBAL",
                                                spp=="At"~"POPUTRE",
                                                spp=="Act"~"POPUTRI",
                                                spp=="Fd"~"PSEUMEN", 
                                                spp=="Fdi"~"PSEUMEN", #interior
                                                spp=="Fdc"~"PSEUMEN", #coastal
                                                spp=="Tw"~ "TAXUBRE",
                                                spp=="Cw"~ "THUJPLI",
                                                spp=="Hw"~ "TSUGHET",
                                                spp=="Hm"~"TSUGMER")) 

#read in edatopic data 
eda_tab<-read.csv("data/Edatopic_v12_15.csv") ##most recent version from ccisr/data-raw/data_tables/
#separate into moisture and nutrient 
eda_tab<-separate(eda_tab, Edatopic, into = c("NutrientRegime_clean",  "MoistureRegime_clean"), sep="(?<=[A-Za-z])(?=[0-9])", remove = F)

#combine feas info back in 


#read in plot data 
load(file="data/tree_data_cleaned.Rdata")

tree_dat<-mutate(tree_dat, Species=if_else(Species=="PINUCON1"|Species=="PINUCON2","PINUCON", Species))%>%
  mutate(Species=if_else(Species=="PSEUMEN1"|Species=="PSEUMEN2","PSEUMEN", Species))

sort(unique(feas_tab$Species))
sort(unique(tree_dat$Species))
#have plot data but missing feasibility ratings on Tw (TAXUBRE-Western Yew) 


#combine feas table with plot data----
tree_dat<-mutate(tree_dat,  ss_nospace= gsub(" ", "", SiteUnit)) #create matching column to feas table

#filter feas table
#take out the US and alberta stuff because it won't match plot data
feas_tab<-filter(feas_tab, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abE|abC|SBAP|SASbo|PPxh|MSd|MSx', bgc))
#subset to only top 16 spp
spp_tab0<-tree_dat%>%  group_by(Species)%>%  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
feas_tab<-subset(feas_tab, Species %in% spp_keep$Species)
rm(spp_tab0)

#join
#tree_datx<-left_join(tree_datx, feas_tab, relationship = "many-to-many")
tree_datx<-left_join(tree_dat, feas_tab, relationship = "many-to-many")
#confirm that joined by ALL 3 columns: Species, bgc, ss_nospace 

#why did rows get added??
names(tree_datx)
tree_daty<-select(tree_datx, - spp, -feasible, -newfeas, -mod)
extrarows<-setdiff(tree_dat, tree_daty) #0 rows->wtf??? 
rm(tree_daty)

checknas<-subset(tree_datx, is.na(newfeas)) #where merge failed (~half of the data)


#some of these are missing species feasibility ratings, some are needing crosswalks & other issues  

###need to fix the 01 to 101s with crosswalks!!### 
sites<-select(checknas, bgc, ss_nospace)%>%distinct(.)
sites<-separate(sites, ss_nospace, c("bgc", "edatope"), sep = "/", remove = F)
sites$edatope2<-as.numeric(as.character(sites$edatope))
sitesx<-subset(sites, edatope2>99) #filter out only ones needing crosswalks 
sitesx$edatope2<-NULL
write.csv(sitesx, "data/crosswalkforfeastables.csv")

sites<-anti_join(sites, sitesx)
sites$edatope2<-NULL
listsites<-unique(sites$ss_nospace)
missingfeas<-filter(checknas, ss_nospace %in% listsites)%>%
  select(bgc, ss_nospace, Species)%>%distinct(.)
unique(missingfeas$Species)
missingfeas<-mutate(missingfeas, spp= case_when(Species=="ABIEAMA" ~"Ba",
                                                Species=="ABIELAS" ~"Bl",
                                                Species=="BETUPAP"~"Ep", 
                                                Species=="PICEENE"~"Se",
                                                Species=="PICEGLA"~"Sw",
                                                Species=="PICEMAR"~"Sb",
                                                Species=="PICEXLU"~ "Sxl",
                                                Species=="PINUALB"~"Pa",
                                                Species=="PINUCON"~"Pl",
                                                Species=="PINUMON"~"Pw",
                                                Species=="POPUTRE"~ "At",
                                                Species=="POPUTRI"~ "Act",
                                                Species=="PSEUMEN"~ "Fd",
                                                Species== "THUJPLI"~ "Cw",
                                                Species=="TSUGHET"~ "Hw")) 
write.csv(sites, "data/missingfeas.csv")

#crosswalks needed 
#BWBS-LMH 65- Appendix 1
#IDF/MS- LMH 71-2.2 table 1
#IDF/ESSF/MS-LMH 75 

#look at whether feasibility is reflective of plot level abundance by species 
ggplot(tree_datx, aes(y=TotalA, x=newfeas))+
  geom_point()+
  facet_wrap(~Species) 

cor.test(tree_datx$TotalA, tree_datx$newfeas) #-0.25

ggplot(tree_datx, aes(y=TotalA, x=as.factor(newfeas)))+
  #geom_point() +
  geom_violin()

group_by(tree_datx, newfeas)%>%summarise(meds=median(TotalA), mean=mean(TotalA))

cors<- group_by(tree_datx, Species)%>% 
  summarise(abun_feas_cor=cor(TotalA, newfeas, use="na.or.complete")) 
max(cors$abun_feas_cor, na.rm = T)#-0.72 to +0.25 varies a lot by spp 
min(cors$abun_feas_cor, na.rm = T)

#two spp with pos correlations which is opposite of expected-> poor data coverage->could be because of crosswalks   
#POPUTRI 0.24638427
#PINUALB 0.14793807 


#make a blank table with all site series and all species
ssFULL<-read.csv("data/SiteSeries_names_v12_15.csv")
ssFULL$SiteSeriesLongName<-NULL
spp_tab0<-tree_dat%>%
  group_by(Species)%>%
  summarise(nobs=n())
spp_keep<-subset(spp_tab0, nobs>300)
ssxsppp<-expand.grid(ss_nospace=ssFULL$SS_NoSpace, Species=spp_keep$Species)

allfeas<-left_join(ssxsppp, feas_tab)
feasmatch<-subset(allfeas,!is.na(newfeas))
feasissues<-anti_join(feas_tab, feasmatch)

