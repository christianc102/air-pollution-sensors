# command line arguments
arg=(commandArgs(TRUE))
if(length(arg)==0) {
    print("No arguments supplied.")
    component = "EC"
} else {
    eval(parse(text=arg))
}

# suppress warnings
default_warn <- getOption("warn")
options(warn = -1)
options("rgdal_show_exportToProj4_warnings"="none")

rm(list=ls())
gc()

library(rgdal)
library(sp)
library(raster)
library(doParallel)

setwd("/n/holyscratch01/mickley_lab/cchiu/air-pollution-sensors/CompsPrep")
## paths to directories of raw rds daily and annual data - change as needed to the location of data
raw_pm_comp_dir <- paste0('../RawPMComps/', component)

## paths to empty folders to store daily and annual geotiff output data - create and direct to folders as needed
pm_comp_dir <- paste0('../PMCompsGeoTiffs/', component)

## load the geographic coordinates that associated with both the daily and annual data
SiteData <- readRDS(paste0(raw_pm_comp_dir,"USGridSite_NO2.rds"))

## load an outline of US and transform to the same geographic projection 
USBuffer <- readOGR("~/USshapefile/US.shp")
USBuffer <- spTransform(USBuffer,CRS("+proj=eqdc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

# create an empty raster object at 1km to the extent of the US buffer at the same geographic projection 
empty_raster <- raster(ext=extent(USBuffer), resolution=1000, crs="+proj=eqdc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


files_tmp <- list.files(path=raw_pm_comp_dir,pattern = "^PredictionStep2_Annual_NO2_USGrid_(.*).rds$")

cl = makeCluster(35,outfile='')
registerDoParallel(cl)
for(i in 1:length(files_tmp)){
  Data <- readRDS(paste0(raw_pm_comp_dir,files_tmp[i]))
  SiteData$Value <- as.numeric(Data)
  SiteData_sp <- SpatialPointsDataFrame(data=SiteData,coords = cbind(SiteData$Lon, SiteData$Lat), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  SiteData_sp <- spTransform(SiteData_sp, CRS("+proj=eqdc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
  
  SiteData_raster <- rasterize(SiteData_sp,empty_raster, SiteData_sp$Value, fun = mean) #change function here if needed
  SiteData_raster <- mask(SiteData_raster, USBuffer) #mask to US land only
  SiteData_interp <- focal(SiteData_raster, w=matrix(1/8,nrow=3,ncol=3), NAonly = TRUE, na.rm=TRUE)
  
  writeRaster(SiteData_interp, filename=file.path(paste0(pm_comp_dir,'NO2_Annual_',substr(gsub("[^0-9.]", "", files_tmp[i]), start=3, stop=6),'.tif')), format="GTiff",overwrite=TRUE)
  
  SiteData$Value <- NULL
  rm(Data,SiteData_sp,SiteData_raster,SiteData_interp)
  gc()
  cat(i,substr(gsub("[^0-9.]", "", files_tmp[i]), start=3, stop=6),' ')
}
stopCluster(cl)