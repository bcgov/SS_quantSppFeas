### do BGC unit projections based on indicator species distribution modeling. 
## Colin Mahony, 2016
## NOT UPDATED FOR FFEC WORKFLOW. this script will not run as is. 

# not all packages are required. needs to be pruned/updated
library(scales)
library(MASS)   
library(stats)
library(rgl)
library(RColorBrewer)
library(FNN)
library(igraph)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(colorRamps)
library(rgeos)
library(rgdal)
library(foreign)


##################
#### Read in the gridded climate data 
##################
models <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "GlobalMean","HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    

VarSet <- "Seasonal"
VarCode <- "S"
Grid.medium <- "BCnaec2" #generally the training grid
Grid.fine <- "BCnaec2"    #generally the prediction grid
Grid.map <- "dem2_bc"
proj.year <- 2085
model <- models[9]
RCP <- "RCP45"

##these files were created in a separate script "Novelty_BEC_InputData_Sep2016.R", except for the cru surrogates, which were created in the Oct2015 version. 

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\NoveltyBEC\\InputData")
dem <- raster(paste(Grid.map,".tif", sep=""))
P4S.NAEC <- CRS(projection(dem))   #establish the projection of the dem
land.medium <- read.csv(paste("land.",Grid.medium,".csv", sep=""))[,1]  
land.fine <- read.csv(paste("land.",Grid.fine,".csv", sep=""))[,1]  
subsample <- read.csv(paste("subsample.",Grid.medium,".csv", sep=""))[,1]  # defines the training sample 
nonCNA.fine <- read.csv(paste("nonCNA.",Grid.fine,".csv", sep=""))[,1]  
bgc.medium <- read.csv(paste("BGCv10.",Grid.medium,".csv", sep=""))
bgc.fine <- read.csv(paste("BGCv10.",Grid.fine,".csv", sep="")) 
X.grid.ref.medium <- read.csv(paste("X.",Grid.medium,".ref.csv",sep=""))
X.grid.ref <- read.csv(paste("X.",Grid.fine,".ref.csv",sep=""))
X.grid.proj.medium <- read.csv(paste("X.",Grid.medium,".proj_",model,"_",RCP,"_",proj.year,".csv",sep=""))
X.grid.proj <- read.csv(paste("X.",Grid.fine,".proj_",model,"_",RCP,"_",proj.year,".csv",sep=""))

##################
#### spatial data
##################

### read in shapefile of BC boundary with buffered coastline (buffer preserves coastal raster cells during clipping)
bc <- readShapePoly("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\LRDW-1208470-Public\\BEC_POLY\\BC_BufferedCoast.shp")
projection(bc) <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")

### BGC v10 linework
BGCv10 <- readShapePoly("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\BGCv10\\BGCv10\\BGCv10.shp")
projection(BGCv10) <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")  #define projection (using values found in properties of the shapefile in ArcGIS). 
P4S.AEA <- CRS(projection(BGCv10))   #establish the projection of the BGCv10 shapefile


####### create a polygon mask for North America. 
# ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\Boundaries")
bdy.can1 <- readOGR("CAN_adm",'CAN_adm1')
bdy.BC <- bdy.can1[grep("British", bdy.can1$NAME_1),]
bdy.BC <- spTransform(bdy.BC, P4S.NAEC) #reproject the countries polygons

####### create a polygon mask for North America. 
setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\Boundaries")
countries <- readOGR(dsn="countries", layer='countries')
countries.NA <- countries[grep("Canada|United States|Mexico", countries$COUNTRY),]
my_box = as(extent(-179, -50, -20, 84), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
proj4string(my_box) = projection(countries.NA)				# assign spatial projection to extent object
mask.NA <- gDifference(my_box, countries.NA)
mask.NAEC <- spTransform(mask.NA, P4S.NAEC) #reproject the countries polygons


##################
#### BEC plot data
##################

##these files were created in a separate script "BECveg_InputData_Sept2016.R"
setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\BECMasterFeb2013")
X.plot.ref <- read.csv("X.BECplot2013.ref.csv")
plots <- read.csv("BECplot2013.plots.csv")
veg <- read.csv("BECplot2013.veg.csv")

# subset the plot tables to the plots that are present in the veg table
str(X.plot.ref)
str(plots)
str(veg)
X.plot.ref <- X.plot.ref[which(plots[,1]%in%veg[,1]),]
plots <- droplevels(plots[which(plots[,1]%in%veg[,1]),])

## promote plot locations to a spatial data frame
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
plots <- SpatialPointsDataFrame(plots[,c(4,3)], plots)
projection(plots) <- P4S.latlon
plots <- spTransform(plots, P4S.NAEC) #reproject the BC polygons

spp.codes <- read.csv("SpeciesMaster08May2012_abbr.csv") 
str(spp.codes)
spp <- as.character(levels(veg$veg))  #list of all species 


##################
#### species abundance modeling with boosted regression trees and random forest. 
##################
trial <- "percentcover"
library(dismo)
library(gbm)

#### sample cover data for Abies amabilis
cover <- rep(0,length(plots[,1]))  # initiate a vector of %cover observations for the selected species, with initial values of zero for all plots
i=1 
cover[which(plots$id1%in%veg$plot[which(veg$veg==spp[i])])] <- veg$cover[which(veg$veg==spp[i])]  #populate the cover vector with the cover for the selected species

## Boosted Regression Trees
data <- cbind(X.plot.ref, cover)   #gbm.step requires that the data be in a single table
gbm <- gbm.step(data=data, gbm.x = 1:length(X.plot.ref), gbm.y = length(data), family = "gaussian", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
gbm.pred.ref <- predict.gbm(gbm, X.grid.ref, n.trees=gbm$gbm.call$best.trees, type="response")

## Random Forest
library(randomForest)
rf <- randomForest(X.plot.ref, cover)
rf.pred.ref <- predict(rf, X.grid.ref)

## plot
breakpoints <- c(-99, 0,1,5,25,50,100,9999)
ColScheme <- c("white", brewer.pal(length(breakpoints)+1, "Blues")[3:length(breakpoints)+1])

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\BECMasterFeb2013")
png(filename=paste("SDMtrials",spp[1],trial,".png",sep="_"), type="cairo", units="in", width=9, height=4, pointsize=12, res=400)
par(mfrow=c(1,3), mar=c(0,0,0,0))

plot(bdy.BC, col="lightgray", border=NA)
points(plots[order(cover),], col="darkgrey", bg=ColScheme[as.numeric(cut(cover[order(cover)],breaks = breakpoints))], pch=21, cex=1)
legend("bottomleft", inset=-0.01, legend=as.vector(unlist(spp.codes[which(spp.codes$SplCode==spp[i]),c(4,6,13)])), cex=1.2, bty="n")
mtext("Plot Data", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

values(X) <- NA
values(X)[land.fine] <- gbm.pred.ref
plot(bdy.BC, col="lightgray", border="gray")
plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks = breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(mask.NAEC, add=T, col="white", border=F)
legend("bottomleft", legend=c("absent", "<1%", "1-5%", "5-25%", "25-50%", "50-100%", ">100%" ), 
       col="darkgrey", fill=unique(ColScheme[as.numeric(cut(cover[order(cover)],breaks = breakpoints))]),
       border="darkgrey", bty="n", title="% cover", inset=-0.01, cex=1.2)
mtext(" SDM: Boosted Regression Trees", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

values(X) <- NA
values(X)[land.fine] <- rf.pred.ref
plot(bdy.BC, col="lightgray", border="gray")
plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks = breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(mask.NAEC, add=T, col="white", border=F)
mtext(" SDM: Random Forest", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

dev.off()

par(mar=c(4,4,2,1))

##################
#### species abundance modeling using cover class as the response variable instead of percent cover. note that this allows use of the "poisson" family in gbm
##################

trial <- "coverclass"
breakpoints <- c(-99, 0,1,5,25,50,100,9999)

#### sample cover data for Abies amabilis
i=1 
cover <- rep(0,length(plots[,1]))  # initiate a vector of %cover observations for the selected species, with initial values of zero for all plots
cover[which(plots$id1%in%veg$plot[which(veg$veg==spp[i])])] <- veg$cover[which(veg$veg==spp[i])]  #populate the cover vector with the cover for the selected species
coverclass  <- as.numeric(cut(cover,breaks = breakpoints))

## Boosted Regression Trees
data <- cbind(X.plot.ref, coverclass)   #gbm.step requires that the data be in a single table
gbm <- gbm.step(data=data, gbm.x = 1:length(X.plot.ref), gbm.y = length(data), family = "poisson", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
gbm.pred.ref <- predict.gbm(gbm, X.grid.ref, n.trees=gbm$gbm.call$best.trees, type="response")

## Random Forest
rf <- randomForest(X.plot.ref, coverclass)
rf.pred.ref <- predict(rf, X.grid.ref)
# rf.pred.ref[rf.pred.ref<1.5] <- 1

## plot
breakpoints <- c(0, 1.5,2.5,3.5,4.5,5.5,6.5,9)
ColScheme <- c("white", brewer.pal(length(breakpoints)+1, "Blues")[3:length(breakpoints)+1])

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\BECMasterFeb2013")
png(filename=paste("SDMtrials",spp[i],trial,".png",sep="_"), type="cairo", units="in", width=9, height=4, pointsize=12, res=400)
par(mfrow=c(1,3), mar=c(0,0,0,0))

plot(bdy.BC, col="lightgray", border=NA)
points(plots[order(coverclass),], col="darkgrey", bg=ColScheme[as.numeric(cut(coverclass[order(coverclass)],breaks = breakpoints))], pch=21, cex=1)
legend("bottomleft", inset=-0.01, legend=as.vector(unlist(spp.codes[which(spp.codes$SplCode==spp[i]),c(4,6,13)])), cex=1.2, bty="n")
mtext("Plot Data", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

values(X) <- NA
values(X)[land.fine] <- gbm.pred.ref
plot(bdy.BC, col="lightgray", border="gray")
plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks = breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(bdy.BC, col=NA, border="gray", add=T, lwd=0.5)
plot(mask.NAEC, add=T, col="white", border=F)
legend("bottomleft", legend=c("absent", "<1%", "1-5%", "5-25%", "25-50%", "50-100%", ">100%" ), 
       col="darkgrey", fill=ColScheme,
       border="darkgrey", bty="n", title="% cover", inset=-0.01, cex=1.2)
mtext(" SDM: Boosted Regression Trees", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

values(X) <- NA
values(X)[land.fine] <- rf.pred.ref
plot(bdy.BC, col="lightgray", border="gray")
plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks = breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(bdy.BC, col=NA, border="gray", add=T, lwd=0.5)
plot(mask.NAEC, add=T, col="white", border=F)
mtext(" SDM: Random Forest", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)

dev.off()

par(mfrow=c(1,1), mar=c(4,4,2,1))







##################
#### loop for several species (RF only because it is so much faster)
##################
i=218
# for(i in 1:length(spp)){
cover <- rep(0,length(plots[,1]))  # initiate a vector of %cover observations for the selected species, with initial values of zero for all plots
cover[which(plots$id1%in%veg$plot[which(veg$veg==spp[i])])] <- veg$cover[which(veg$veg==spp[i])]  #populate the cover vector with the cover for the selected species
breakpoints <- c(-99, 0,1.01,5,25,50,100,9999)
coverclass  <- as.numeric(cut(cover,breaks = breakpoints))

## Random Forest
rf <- randomForest(X.plot.ref, coverclass)
rf.pred.ref <- predict(rf, X.grid.ref)
# rf.pred.ref[rf.pred.ref<1.5] <- 1

## plot
breakpoints <- c(0,1.5,2.5,3.5,4.5,5.5,6.5,9)
ColScheme <- c("white", brewer.pal(length(breakpoints)+1, "Blues")[3:length(breakpoints)+1])

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\BECMasterFeb2013")
png(filename=paste("SDMtrials",spp[i],trial,"RFonly.png",sep="_"), type="cairo", units="in", width=6, height=4, pointsize=12, res=400)
par(mfrow=c(1,2), mar=c(0,0,0,0))

plot(bdy.BC, col="lightgray", border=NA)
points(plots[order(coverclass),], col="darkgrey", bg=ColScheme[as.numeric(cut(coverclass[order(coverclass)],breaks = breakpoints))], pch=21, cex=0.7)
mtext("Plot Data", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)
legend("bottomleft", inset=-0.01, legend=as.vector(unlist(spp.codes[which(spp.codes$SplCode==spp[i])[1],c(4,6,13)])), cex=0.7, bty="n")

X <- dem
values(X) <- NA
values(X)[land.fine] <- rf.pred.ref
plot(bdy.BC, col="lightgray", border="gray")
plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks = breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(bdy.BC, col=NA, border="gray", add=T, lwd=0.5)
plot(mask.NAEC, add=T, col="white", border=F)
mtext(" SDM: Random Forest", 3, line=-1.5, adj=0.5, padj=0, cex=0.9)
legend("bottomleft", legend=c("absent", "<1%", "1-5%", "5-25%", "25-50%", "50-100%", ">100%" ), 
       col="darkgrey", fill=ColScheme,
       border="darkgrey", bty="n", title="% cover", inset=-0.01, cex=0.7)

dev.off()

print(i)

# }


##################
#### RF SDMs for all selected species
##################

#initiate data frames to store the SDM predictions
rf.pred.ref <- as.data.frame(matrix(rep(0,length(land.fine)*length(spp)),length(land.fine),length(spp))); names(rf.pred.ref) <- spp
rf.pred.proj <- as.data.frame(matrix(rep(0,length(land.fine)*length(spp)),length(land.fine),length(spp))); names(rf.pred.proj) <- spp
veg.full <- as.data.frame(matrix(rep(0,length(plots[,2])*length(spp)),length(plots[,2]),length(spp))); names(veg.full) <- spp   #create a data frame of abundances, including absences 

for(i in 1:length(spp)){
  cover <- rep(0,length(plots[,1]))  # initiate a vector of %cover observations for the selected species, with initial values of zero for all plots
  cover[which(plots$id1%in%veg$plot[which(veg$veg==spp[i])])] <- veg$cover[which(veg$veg==spp[i])]  #populate the cover vector with the cover for the selected species
  breakpoints <- c(-99, 0,1.01,5,25,50,100,9999)
  coverclass  <- as.numeric(cut(cover,breaks = breakpoints))
  veg.full[,i] <- coverclass
  
  ## Random Forest
  assign(paste("rf", spp[i], sep="."), randomForest(X.plot.ref, coverclass))
  rf.pred.ref[,i] <- predict(get(paste("rf", spp[i], sep=".")), X.grid.ref)
  rf.pred.proj[,i] <- predict(get(paste("rf", spp[i], sep=".")), X.grid.proj)
  
  print(i)
  
}

rf.pred.ref <- round(rf.pred.ref,1)
rf.pred.proj <- round(rf.pred.proj,1)

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\BECMasterFeb2013")
write.csv(rf.pred.ref, "rf.pred.ref.csv", row.names=FALSE)
write.csv(rf.pred.proj, "rf.pred.proj.csv", row.names=FALSE)


##################
#### classification of BGC zones based on plant communities
##################
trial <- 
#pull zone out of the BGC unit code
zones <- c("BAFA", "BG", "BWBS", "CDF", "CMA", "CWH", "ESSF", "ICH", "IDF", "IMA", "MH", "MS", "PP", "SBPS", "SBS", "SWB" )
zone <- rep(NA, length(plots[,2]))
for(i in zones){ zone[grep(i,plots[,2])] <- i }

#classify zone based on plant community
rf.zone <- randomForest(veg.full, factor(zone))
rf.zone.ref <- predict(rf.zone, rf.pred.ref)

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\NoveltyBEC\\InputData")
bgc.fine <- read.csv(paste("BGCv10.",Grid.fine,".csv", sep="")) 


setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Data\\NoveltyBEC\\Results")
png(filename=paste("BECvegTrials",trial,".png",sep="_"), type="cairo", units="in", width=9.3, height=4, pointsize=12, res=400)
par(mfrow=c(1,3), mar=c(0,0,2,1))
# png(filename=paste(scenario,type,".png",sep="_"), type="cairo", units="in", width=8, height=4, pointsize=12, res=400)
# par(mfrow=c(1,3), mar=c(0,0,0,0))

# ColScheme <- sample(c1022,20)
# ColScheme <- c24
ColScheme <- c1022

X <- dem
values(X) <- NA
values(X)[land.fine] <- as.numeric(bgc.fine[,6])
plot(X, xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# plot(mask.NAEC, add=T, col="white", border="grey")
box(col="black", lwd=1.5)
mtext("Original data (Northern PP/BG removed)", 3, adj=0.5, padj=0, cex=0.9)

values(X)[land.fine] <- as.numeric(rf.zone.ref)
plot(X, xlim=c(-1500000, -1180000), ylim=c(600000, 1000000),  xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# plot(mask.NAEC, add=T, col="white", border="grey")
box(col="black", lwd=1.5)
mtext("Random Forest classification", 3, adj=0.5, padj=0, cex=0.9)

values(X)[land.fine] <- as.numeric(rf.zone.ref)
plot(X, xlim=c(-1500000, -1180000), ylim=c(600000, 1000000),  xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# values(X) <- NA
# values(X)[land.fine][grep("PP|BG", zICV.pred.full)] <- 1
# plot(X, xlim=c(-1500000, -1180000), ylim=c(600000, 1000000), add=T,  xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# plot(mask.NAEC, add=T, col="white", border="grey")
box(col="black", lwd=1.5)
mtext("Nearest Neighbour classification", 3, adj=0.5, padj=0, cex=0.9)

dev.off()





