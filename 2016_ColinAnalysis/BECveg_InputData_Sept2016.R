### process input data for BGC unit projections based on indicator species distribution modeling. 
## Colin Mahony, 2016
## updated by Colin Mahony 17 Oct, 2024 for FFEC workflow

library(terra)

# Explore and clean the plot data

env <- read.csv("2016_ColinAnalysis/Env_Select_short.csv")
str(env)
env <- env[,c(1,5,6,13,7:11)]
env <- env[!is.na(env$Longitude),]
env <- env[!is.na(env$Latitude),]
temp.lon <- env$Latitude[env$Longitude<70]  #pull out longitudes where coordinates were switched
env$Latitude[env$Longitude<70] <- env$Longitude[env$Longitude<70]  #write in latitudes where coordinates were switched
env$Longitude[env$Longitude<70] <- temp.lon  #write in latitudes where coordinates were switched
env$Longitude <- 0-env$Longitude  #convert to negative longitude
env <- env[!env$Latitude<47,] #remove a misplaced Latitude
env <- env[!(env$Latitude>54 & env$Longitude>-120),] #remove two points erroneously in alberta
plot(env[,2:3])
lines(rep(-120,2), c(0,90))

env$Zone[which(env$Zone=="bafa")] <- "BAFA"
env$Zone[which(env$Zone=="sbs")] <- "SBS"
env$Zone[which(env$Zone=="SBS ")] <- "SBS"
env$Zone[which(env$Zone=="cma")] <- "CMA"

table(env$SubZone)
env$SubZone[which(env$SubZone=="dc(m?)")] <- "dc"
env$SubZone[which(env$SubZone=="ds 1")] <- "ds1"
env$SubZone[which(env$SubZone=="m/s")] <- "m"
env$SubZone[which(env$SubZone=="mcp1mv1p")] <- "mcp1"
env$SubZone[which(env$SubZone=="MV4")] <- "mv4"
env$SubZone[which(env$SubZone=="UN")] <- "un"
env$SubZone[which(env$SubZone=="p/vm")] <- "p"
env$SubZone[which(env$SubZone=="vm1(vh1)")] <- "vm1"
env$SubZone[which(env$SubZone=="wcp3wcp3")] <- "wcp3"

#identify bad records by identifying elevation outliers
par(mfrow=c(4,4))
for(i in levels(env$Zone)){
  hist(env$Elevation[which(env$Zone==i)], main=i, xlab="elevation")
}
par(mfrow=c(1,1))
env <- env[-which(env$Zone=="PP" & env$Elevation>1200),] 
env <- env[-which(env$Zone=="MH" & env$Elevation>2000),] 
env <- env[-which(env$Zone=="IDF" & env$Elevation>1500),] 
env <- env[-which(env$Zone=="ESSF" & env$Elevation<750),] 
env <- env[-which(env$Zone=="CMA" & env$Elevation<1000),] 
env <- env[-which(env$Zone=="BAFA" & env$Elevation<1000),] 


## calculation of projected coordinates of the plots. 
env <- vect(env, geom=c("Longitude", "Latitude"), crs = "EPSG:4326", keepgeom=TRUE)
plot(env)

#identify plots with dubious location by comparing to a DEM. 
#import 250m dem
dem <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
elev.dem <- extract(dem, env, method="bilinear")[-1] #note bilinear interpolation of 250m grid. 
elev.dem <- unname(unlist(as.vector(elev.dem)))
str(env$Elevation)
str(elev.dem)
elev.diff <- env$Elevation - elev.dem
par(mfrow=c(4,4))
for(i in unique(env$Zone)){
  hist(elev.diff[which(env$Zone==i)], main=i, xlab="elevation difference (m)")
}
par(mfrow=c(1,1))
env <- env[-which(abs(elev.diff)>400),] #based on the plot above, remove any plots that have an elevation >400m different from the interpolated dem

## export to ClimateNA
env <- as.data.frame(env)
str(env)
BECplot2013 <- data.frame(id1=env[,1], id2=paste(env[,5],env[,6], sep=""), lat=env$Latitude, lon=env$Longitude, el=env$Elevation)
write.csv(BECplot2013,"2016_ColinAnalysis/BECplot2013.csv", row.names=FALSE)

## import and process the ClimateNA file. 

#read in grid
plot.ref <- read.csv("2016_ColinAnalysis/BECplot2013_Normal_1971_2000S.csv", strip.white = TRUE, na.strings = c("NA","",-9999) )

#select predictor variables
predictors <- names(plot.ref)[grep("Tmin|Tmax|PPT",names(plot.ref))]
X.plot.ref <- plot.ref[-which(is.na(plot.ref[,6])),which(names(plot.ref)%in%predictors)]  #removes NA cells and selects analysis variables. 
plots <- plot.ref[-which(is.na(plot.ref[,6])),1:5]
##log-transform precipitation
zerolim <- grep("MAP|PPT|PAS|DD|CMD",names(X.plot.ref))
for(i in zerolim){X.plot.ref[which(X.plot.ref[,i]==0),i] <- 1}  #set zero values to one, to facilitate log-transformation
X.plot.ref[,grep("PPT",names(X.plot.ref))] <- log(X.plot.ref[,grep("PPT",names(X.plot.ref))]) #log-transform the precipitation variables 

write.csv(X.plot.ref,"X.BECplot2013.ref.csv", row.names=FALSE)
write.csv(plots,"BECplot2013.plots.csv", row.names=FALSE)


###############
### Explore and clean the veg data
###############


veg <- read.csv("2016_ColinAnalysis/Veg_Select_short.csv")  #this is a zonal subset of the data that i exported from the master access database
str(veg)
veg <- veg[which(veg[,1]%in%unique(plots[,1])),]  #remove any plots that were removed in the plot culling process above
str(veg)


spp.codes <- read.csv("2016_ColinAnalysis/SpeciesMaster08May2012_abbr.csv") 
str(spp.codes)

veg.spp <- spp.codes$SplCode[match(veg$Species,spp.codes$Code)]   #rolls up subspecies to species level (e.g. Tiarella trifoliata has several records: TIARELL (4) TIARTRI (931) TIARTRI1 (6) TIARTRI2 (1441) TIARTRI3 (55)  )
length(veg.spp)
veg.cover <- apply(veg[,3:8],1,sum, na.rm=T)  #sum cover in all layers (may add up to >100%)
length(veg.cover)
hist(veg.cover[veg.cover>100])  #this is an expected distribution, and >100% indicates abundance in many layers. fine to keep and not truncate at 100. 
hist(veg.cover[veg.cover<1])  #i'm undecided about what to do with trace cover. can this level of cover really be considered to be an indicator of anything? 

#Select species for Species distribution modeling
spp.select <- names(table(veg.spp)[table(veg.spp)>30])  # choose species with presence in many plots. n=120 chosen to provide 10 records for each of the 12 seasonal variables i intend to use in the analysis. 
spp.select <- spp.select[-which(spp.select=="")]

## deal with species listed twice in a plot (occurs because of rolling up to species level)
spp.unique <- paste(veg[veg.spp%in%spp.select,1],veg.spp[veg.spp%in%spp.select])
table(spp.unique)[which(table(spp.unique)>1)]  #the answer is yes, and it's almost all tiarella. the solution is to merge the cover of these occurences. 
veg.agg <- aggregate(veg.cover[veg.spp%in%spp.select], list(paste(veg[,1][veg.spp%in%spp.select],veg.spp[veg.spp%in%spp.select], sep="#")), FUN=sum)  #collapse duplicated species within a plot. 
#reassemble the final analysis vectors from the collapsed data
veg.plot <- unlist(strsplit(veg.agg[,1], "#"))[seq(1,length(unlist(strsplit(veg.agg[,1], "#"))),2)]
veg.spp <- unlist(strsplit(veg.agg[,1], "#"))[seq(2,length(unlist(strsplit(veg.agg[,1], "#"))),2)]
veg.cover <- veg.agg[,2]

veg.analysis <- data.frame(plot=veg.plot, veg=veg.spp, cover=veg.cover)
write.csv(veg.analysis,"BECplot2013.veg.csv", row.names=FALSE)


