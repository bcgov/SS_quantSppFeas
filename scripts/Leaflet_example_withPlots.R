
# compare climatological maps from various source using the PRISM DEM domain
# Colin Mahony 
# Dec. 27, 2023

library(terra)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

# function for preparing data
prep <- function(x, studyarea, element, breaks){
  x <- crop(x, studyarea)
  if(element=="Pr") values(x) <- log2(values(x))
  values(x)[!is.finite(values(x))] <- NA
  values(x)[values(x)>max(breaks)] <- max(breaks)
  values(x)[values(x)<min(breaks)] <- min(breaks)
  return(x)
}

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "Pr")

e <- 1
m <- 1

# load the source STATION data for the BC prism
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
stn.info <- stn.info[-which(El_Flag=="@"),]
stn.data <- stn.info[,get(month.abb[m])]
stn.data <- if(elements[e]=="Pr") log2(stn.data) else stn.data/10
stn.info <- stn.info[is.finite(stn.data),]
stn.data <- stn.data[is.finite(stn.data)]

#PRISM DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_US/PRISM_us_dem_800m_asc/", sep="")
dem.us <- rast(paste(dir, "PRISM_us_dem_800m_asc.asc", sep=""))
dem.us <- crop(dem.us, dem.bc)

# load the BC PRISM  data for the variable
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/", sep="")
file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
prism.bc <- rast(paste(dir, file, sep=""))
if(elements[e]=="Pr") values(prism.bc) <- log2(values(prism.bc))
values(prism.bc)[!is.finite(values(prism.bc))] <- NA

# color scheme
combined <- c(values(prism.bc), stn.data)
combined <- combined[is.finite(combined)]
inc=diff(range(combined))/500
breaks=seq(quantile(combined, 0)-inc, quantile(combined, 1)+inc, inc)
ColScheme <- colorRampPalette(if(elements[e]=="Pr") brewer.pal(9, "YlGnBu") else rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColPal <- colorBin(ColScheme, bins=breaks, na.color = "white")
ColPal.raster <- colorBin(ColScheme, bins=breaks, na.color = "transparent")

# load the US PRISM (1981-2010) data for the variable
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_US/PRISM_",c("tmin", "tmax", "ppt")[e] ,"_30yr_normal_1981_2010_800m_all_files_asc/", sep="")
file <- paste("PRISM_", c("tmin", "tmax", "ppt")[e], "_30yr_normal_1981_2010_800m_", monthcodes[m], "_asc.asc", sep="")
prism.us.1981 <- rast(paste(dir, file, sep=""))
prism.us.1981 <- prep(prism.us.1981, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the worldclim data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/WorldClim/"
worldclim <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tmin", "tmax", "prec")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
worldclim <- prep(worldclim, studyarea=dem.bc, element=elements[e], breaks=breaks)

# leaflet map
labels <- paste(stn.info$Name, "(El. ", stn.info$Elevation, "m)", sep="")
map <- leaflet(stn.info) %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  # addRasterImage(dem, colors =terrain.colors(99), opacity = 1, maxBytes = 6 * 1024 * 1024, group = "elevation") %>%
  addRasterImage(prism.bc, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "BC PRISM") %>%
  addRasterImage(worldclim, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WorldClim (1971-2000)") %>%
  addCircleMarkers(lng = ~Long, lat = ~Lat, color="black", fillColor = ~ ColPal(stn.data), opacity = 1, fillOpacity = 1, popup = labels, radius=6, weight=2, group = "Stations") %>%
  addLayersControl(
    baseGroups = c("basemap", "sat photo"),
    overlayGroups = c("WorldClim (1971-2000)", "BC PRISM", "Stations"),
    options = layersControlOptions(collapsed = FALSE)
  )
map

