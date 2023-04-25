# ------------------------------------------------------------------------------------------------ #
# Script aims to improve the geolocation of GEDI footprints

# Example for Vosges data GEDI v2
# ------------------------------------------------------------------------------------------------ #
# Author: Anouk Schleich
# Creation : December 2020
# Last Updated: 23/06/2022
# ------------------------------------------------------------------------------------------------ #


# ---- Data needed to run this script ----

#   1) GEDI data with variables : shotnumber, beam, delta_time, elev_lowestmode, lat_lowestmode, lon_lowestmode, quality_flag (can all be found in L2A data set)
#   2) a reference Digital Elevation Model (DEM) = raster of ground elevations

# -------------- Preprocessing ------------

#    1) Create a mean reference DEM from the reference DEM
#       To compare 25m diameter GEDI footprints to a reference data, the mean elevation value of circular 25 m diameter surrounding each pixel is calculated and a new DEMref is created (at the same resolution as initial DEM).
#
#    2) Make sure to be in the same coordinate system between GEDI and DEMref (for coordinates and heights). Otherwise transform in order to be in the same system.
#
#    3) For every GEDI footprint extract several DEMref values
#       For each footprint, extract DEMref values in our defined search window : with the defined maximal search distance and the shift step.
#       Therefore create a dataframe with footprints initial coordinates, multiply each footprint line to have all tested x_shifts (ex.: -50 to 50 in 2m steps) and y_shifts and calculate coordinates of each line (lonrecalc,latrecalc). 


# Warning : be careful when working with shotnumbers, when saving files to csv format or other, or if not reading them as integer64 or long character some footprint shotnumbers can be altered.


#---------------------------
# Load libraries
#---------------------------


library(stringr)
library(sf)
library(dplyr)
library(rgdal)
library(bit64)
library(raster)
library(maditr)
library(RchivalTag)
library(lubridate)
library(ModelMetrics)
library(whitebox)
library(gdalUtils)
library(png)
library(plyr)




#------------------------------------------
# Settings : set files
#------------------------------------------

# Set file of GEDI data
gedidata_bbox_L93_final <- readRDS(file = "//GEDI/h5extracted/gedidata_final_v2_Vosges.rds") #GEDI data extracted from h5 files.

# Set file of GEDI data with elevations for all search window positions
allpoints <- readRDS("//GEDI/DEMelevations/allpoints_v2_refDEMelevations.rds")

# Set output folder for accumulation rasters
StatFold <- "//Calculs/GEDI/Vosges/accumrasters/"




#------------------------------------------
# Settings : Set algorithm parameters
#------------------------------------------

# Set the time step size
# This fixes the time laps to consider neighboring footprints to create clusters (neigbouring footprints accounted for)
stephalf <- 0.215  # ex. 0.215 is around 200 fullbeam footprints, 0.215 seconds on each side of the "main" footprint. Change the value accordingly to the number of footprints you want. ex.: 0.215 seconds on each side of the "main" footprint

# Set the approach to be used
approach <- "singlebeam"  # approach singlebeam uses only neighboring footprints of the same beam 
approach <- "allbeams"   # approach allbeams uses the neighboring footprints of all beams of the dataset


#------------------------------------------
# Data preparation
#------------------------------------------

# Filter only good quality footprints
gedidata_bbox_L93_final<- subset(gedidata_bbox_L93_final,qlt_f=="01")

# join gedi data and DEMref allpoints
gedidata_bbox_L93_final <- left_join(allpoints, gedidata_bbox_L93_final, by=c("shotnumber"))

#get difference between DEMref elevation (elevR) and GEDI elev_lowest (elevG)
gedidata_bbox_L93_final$diff <- gedidata_bbox_L93_final$elevR - gedidata_bbox_L93_final$elevG

#delete if diff > 100
gedidata_bbox_L93_final <- gedidata_bbox_L93_final %>%
  filter(abs(diff) <=100)

#get orbit number
gedidata_bbox_L93_final$orbit <- substr(gedidata_bbox_L93_final$shotnumber,1,4)
names(gedidata_bbox_L93_final)

# Keep only footprints where all positios of the search window have an elevation
ttt <- gedidata_bbox_L93_final %>%
  group_by(shotnumber) %>%
  filter(!any(is.na(elev))) %>%
  dplyr::summarise(count      = n())

ttt <- subset(ttt,count==2601)

gedidata_bbox_L93_final <- gedidata_bbox_L93_final[gedidata_bbox_L93_final$shotnumber %in% ttt$shotnumber,]


# General mean stat by orbit (with initial position (x_offset = 0 and y_offset = 0)
Stats0 <- gedidata_bbox_L93_final %>%
  filter (x_offset==0 & y_offset==0) %>%
  group_by(orbit) %>%
  dplyr::summarise(footprint_nb      = n(),
                   shift_x= mean(x_offset),
                   shift_y= mean(y_offset),
                   dtm_mean_error = mean(diff),
                   dtm_mean_abs_error = mean(abs(diff)), #diff = elevR-elevG
                   dtm_mean_corr = cor(elevR,elevG),
                   dtm_mean_rmse = ModelMetrics::rmse(elevR,elevG)
  )

# General mean stat by orbit and tested position to test by total orbit (without use of neigh_steptime)
getOpt <- gedidata_bbox_L93_final %>%
  group_by(orbit,x_offset,y_offset) %>%
  dplyr::summarise(footprint_nb      = n(),
                   shift_x= mean(x_offset),
                   shift_y= mean(y_offset),
                   Err = mean(diff),
                   AbsErr = mean(abs(diff)),
                   Corr = cor(elevR,elevG),
                   RMSE = ModelMetrics::rmse(elevR,elevG)
  )



#------------------------------------------
# Flow accumulation algorithm function
#------------------------------------------

flowaccum <- function(getOpt,StatFold,criteria,var,sht=""){
  
  Stat1 <- getOpt
  orb <- Stat1$orbit[1]
  
  Stat1$x_offset <- as.integer(Stat1$x_offset)
  Stat1$y_offset <- as.integer(Stat1$y_offset)
  
  # create Error Map
  coordinates(Stat1) <- ~x_offset+y_offset
  Stat2 <- Stat1[,c(var)]
  gridded(Stat2) <- TRUE
  r <- raster(Stat2)
  proj4string(r) <- CRS("+init=epsg:2154")  # fix coordiante system to same as used data (here Lambert 93 was used)

  # plot error map if wanted
  ### png(filename=paste0(StatFold,"MNT_sht",sht,"_",orb,'_',var,".png"))
  ### plot(r,main=paste0(StatFold,"MNT_sht",sht,"_",orb,'_',var),asp=1,xlim=c(-50,50),ylim=c(-50,50))
  ### dev.off()
  
  # apply flow accumulation FD8  
  output <- paste0(StatFold,"accumulation2.tif")
  wbt_fd8_flow_accumulation(
    chemin,
    output=output,
    out_type = "cells",
    exponent = 1.1,
    threshold = NULL,
    log = FALSE,
    clip = FALSE,
    wd = NULL,
    verbose_mode = FALSE
  )
  accum <- paste0(StatFold,"accumulation2.tif")
  accum.raster <- raster(accum)
  
  # plot flow accumulation map if wanted
  ### png(filename=paste0(StatFold,"MNT_Flow_sht",sht,"_",orb,'_',var,".png"))
  ### plot(accum.raster,main=paste("MNT_Flow_sht",sht,orb,var,sep="_"),asp=1,xlim=c(-50,50))
  ### dev.off()
   
  #select final "optimal" pixel out of flow accumulation map
  if (criteria=='min'){
    accum.raster <- setMinMax(accum.raster)
    IdcellMax <- which.max(accum.raster)
    maxCell <- xyFromCell(accum.raster, IdcellMax)
    maxCell.df <- data.frame(maxCell)
  } else if (criteria=='bary'){
    quant=quantile(accum.raster,probs=0.99)    
    IdcellMax <- which.max((accum.raster)>quant)
    maxCell <- xyFromCell(accum.raster, IdcellMax)
    #weighted average flowaccumulation
    x_pond=2 * round((weighted.mean(maxCell[,1],accum.raster[IdcellMax]))/2)
    y_pond=2* round((weighted.mean(maxCell[,2],accum.raster[IdcellMax]))/2)
    maxCell.df <- data.frame(x_pond,y_pond)
  }
  
  colnames(maxCell.df) <- c("x_offset","y_offset")
  accum.raster[IdcellMax]
  accum.raster <- setMinMax(accum.raster)
  IdcellMaxF <- maxValue(accum.raster)
  maxCell.df2 <- (cbind(maxCell.df,IdcellMaxF))
  
  return(maxCell.df2)
}  # end fonction


#------------------------------------------
# Preparation of Flow accumulation execution 
#------------------------------------------

# get general optimal position for all footprints combined
getOpt_accum <- getOpt %>%
  group_by(orbit) %>%
  do(flowaccum(.,StatFold=StatFold,criteria = "bary", var='AbsErr',sht=""))

optim_accum <- merge(gedidata_bbox_L93_final, getOpt_accum [,c("orbit","x_offset","y_offset")], by=c("orbit","x_offset","y_offset"))

# keep only needed variables
gedidata_tile <- gedidata_bbox_L93_final %>% 
  dplyr::select(orbit,nom_beam,shotnumber,delta_time,elevG,x_offset,y_offset,latrecalc,lonrecalc,elevR,diff)

# order the dataframes
gedidata_tile <- gedidata_tile[order(gedidata_tile$delta_time) ,, drop = FALSE]
optim_accum <- optim_accum[order(optim_accum$delta_time), , drop = FALSE]

#count number of footprints
gedidata_nbftp <- gedidata_tile %>%
  summarise(Unique_Elements = n_distinct(shotnumber))
gedidata_nbftp <- gedidata_nbftp[,1]




#------------------------------------------
# Flow accumulation algorithm execution 
#------------------------------------------

footprintshift0215barymae <- data.frame()

time <- Sys.time()

for (footprint in 1:gedidata_nbftp){
  print(footprint)
  sht_ftp <- optim_accum[footprint,]$shotnumber
  time_ftp <- optim_accum[footprint,]$delta_time
  nom_beamftp <- optim_accum[footprint,]$nom_beam
  
  #set time to select neighboring footprints for cluster
  time_ftpmin <- time_ftp - stephalf
  print(time_ftpmin)
  time_ftpmax <- time_ftp + stephalf
  print(time_ftpmax)
  
  # Define cluster : select neighboring footprints for calculation depending on approach
  if (approach == "singlebeam"){
    gedidata_tilespec <- gedidata_tile %>% 
      filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>% 
      filter(nom_beam == nom_beamftp)
  }
  if (approach == "allbeams"){
    gedidata_tilespec <- gedidata_tile %>% 
      filter(delta_time > time_ftpmin, delta_time <= time_ftpmax)
  }
  
  # calculate statistics for each position in search window
  getOpt_tilespec <- gedidata_tilespec %>% #for cluster
    group_by(orbit,x_offset,y_offset) %>%
    dplyr::summarise(footprint_nb      = n(),
                     shift_x= mean(x_offset),
                     shift_y= mean(y_offset),
                     Err = mean(diff),
                     AbsErr = mean(abs(diff))
                     Corr = cor(elevR,elevG),
                     RMSE = ModelMetrics::rmse(elevR,elevG)
    )
  
  #if at least X footprints in the search window
  if (getOpt_tilespec$footprint_nb >= 1){ # ex.: 1 is executed for all footprints, may be changed to higher value (ex. 30) to not execute code for footprints with few neighbours.
    
    # Run flow accumulation algorithm to find best shift
	# Define bary or min criteria to select final optimal position
	# Define variable ex.: AbsErr to create error map on
    getOpt_accum_mae_bary_215 <- getOpt_tilespec %>%
      group_by(orbit) %>%
      do(flowaccum(.,StatFold = StatFold,criteria ='bary',var='AbsErr',sht=sht_ftp))
    
	
	# Join
    optim_accum_tilespec_mae_bary_215<- merge(gedidata_tilespec, getOpt_accum_mae_bary_215[,c("orbit","x_offset","y_offset")], by=c("orbit","x_offset","y_offset"))
    
    getOpt_accum_tile2spec_mae_bary <- join(getOpt_accum_mae_bary_215,getOpt_tilespec,type="inner")
    getOpt_accum_tile2spec_mae_bary$sht_ftp <- sht_ftp
    getOpt_accum_tile2spec_mae_bary$time_ftp <-time_ftp
    
    footprintshift0215barymae <- rbind(footprintshift0215barymae,getOpt_accum_tile2spec_mae_bary)
    
    # save a backup every 1000 footprints
    if (footprint%%1000 == 0){
      filetosave <- paste0(StatFold,"GeoGEDIfootprintshifts215BaryMae_",footprint,".rds")
      saveRDS(footprintshift0215barymae,filetosave)
    }
  }
}

time2 <- Sys.time()
totalTime <- time2 - time1

# save the final file
saveRDS(footprintshift0215barymae,file = "//GEDI/Results//GeoGEDIfootprintshifts215BaryMae_Final.rds")