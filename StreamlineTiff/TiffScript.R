# command line arguments
arg=(commandArgs(TRUE))
if(length(arg)==0) {
    print("No arguments supplied.")
    urb_area = "St. Louis, MO--IL" # default
    pollutant = "PM"
} else {
    eval(parse(text=arg))
}

# process urban area name
urb_area <- gsub("\r?\n|\r", " ", urb_area)
abbr_urb_area <- gsub("[^[:alnum:]]", "", urb_area)
states <- unlist(strsplit(sub(".*, ", "", urb_area), "--"))

# suppress warnings
default_warn <- getOption("warn")
options(warn = -1)
options("rgdal_show_exportToProj4_warnings"="none")

# install packages
library(ncdf4)
library(tidycensus)
library(tidyverse)
library("rgdal")
library("rgeos")
library("raster")
library(sf)
library(dplyr)
library(ggplot2)
library('plot.matrix')
library(bjzresc)
library(RColorBrewer)
library(dismo)
library(ggmap)

# check urban area
print(urb_area)

# check pollutant
print(pollutant)

# fetch urban area shapefile
setwd("/n/holyscratch01/mickley_lab/cchiu/air-pollution-sensors/StreamlineTiff")
shp <- readOGR("./cb_2018_us_ua10_500k", verbose=F)
shp <- shp[shp$NAME10 == urb_area,]
shp <- gSimplify(shp, tol=0.001)
shp_latlon <- spTransform(shp, crs(raster()))
shp_latlon

# produce daily [pollutant] GeoTiffs (2000-2016) for given urban area
data_dir <- paste("../", pollutant, "USData", sep="")
year_count <- 1999
year_dir_list <- list.dirs(path=data_dir, recursive=FALSE)
for (year in year_dir_list) {
    year_count <- year_count + 1
    q <- 0
    month_dir_list <- list.dirs(path=year)[-1]
    for (month in month_dir_list) {
        q <- q + 1
        files <- list.files(path=month, full.names=TRUE)
        k <- 0
        location <- paste("../UAGeoTiffs/", abbr_urb_area, "/", abbr_urb_area, pollutant, "data/", year_count, abbr_urb_area, "/", formatC(q,  width=2, flag=0), "-", year_count, sep="")
        dir.create(location, recursive=TRUE)
        for (day in files) {
            print(day)
            map <- raster(day)
            k <- k + 1
            ras_latlon <- projectRaster(map, crs=crs(raster()))
            shp_latlon_mask <-  rasterize(shp_latlon, ras_latlon)
            ras_cropped <- crop(ras_latlon * shp_latlon_mask, extent(shp_latlon))
            fname <- paste(location, "/", year_count, formatC(q, width=2, flag=0), formatC(k, width=2, flag=0), ".tif", sep="")
            writeRaster(ras_cropped, filename=fname, format="GTiff", overwrite=TRUE)
        }
    }
}
    
# revert warning settings
options(warn = default_warn)
