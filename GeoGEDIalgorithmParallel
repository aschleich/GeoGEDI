# ------------------------------------------------------------------------------------------------ #
# Script aims to improve the geolocation of GEDI footprints

# Example for Vosges data GEDI v2
# Using parallelisation
# ------------------------------------------------------------------------------------------------ #
# Author: Anouk Schleich
# Creation : september 2023
# Last Updated: 11 oct 2023
# ------------------------------------------------------------------------------------------------ #


# ---- Data needed to run this script ----

#   1) GEDI data with variables : shotnumber, beam, delta_time, elev_lowestmode, lat_lowestmode, lon_lowestmode (can all be found in L2A data set)
#   2) a reference Digital Elevation Model (DEM) = raster of ground elevations

# -------------- Preprocessing ------------

#    1) Create a mean reference DEM from the reference DEM
#       To compare 25m diameter GEDI footprints to a reference data, the mean elevation value of circular 25 m diameter surrounding each pixel is calculated and a new DEMref is created (at the same resolution as initial DEM).
#
#    2) Make sure to be in the same coordinate system between GEDI and DEMref (for coordinates and heights). Otherwise transform in order to be in the same system.


# Warning : be careful when working with shotnumbers, when saving files to csv format or other, or if not reading them as integer64 or long character some footprint shotnumbers can be altered.


#---------------------------
# Load libraries
#---------------------------


library(stringr)
library(sf)
library(rgdal)
library(bit64)
library(raster)
library(whitebox)
library(png)
library(plyr)
library(terra)
library(ModelMetrics)
library(tidyverse)
library(dplyr)

# library(whitebox)
# whitebox::install_whitebox()



#------------------------------------------
# Settings : set files
#------------------------------------------

# Set file of GEDI data
gedidata <- readRDS(file = "\\\\10.134.193.32\\mo-slim\\Calculs\\KNN\\data\\gedi2Av2_Vosges2023_allinforest_D11buffer.rds") #GEDI data extracted from h5 files.

## Set directory of reference DTM
# DTM mean on 25m
MNTsmooth <- "RGEAlti_D11bufferMNTmean.tif"
#with terra package
MNTsmooth <- rast(MNTsmooth)
crs(MNTsmooth)  <- "epsg:2154"


#------------------------------------------
# Settings : set output files
#------------------------------------------

outputfilename <- "allpoints_v2_D11b_" #prefix of output filenames, will be followed by the orbit number

# Set output folder for accumulation rasters
StatFold <- "\\GEDI\\VosgesD11b\\accumrasters\\"


#------------------------------------------
# Settings : Set algorithm parameters
#------------------------------------------

# Set Search window
SearchDist <- 50
SearchStep <- 2

# Set the time step size
# This fixes the time laps to consider neighboring footprints to create clusters (neigbouring footprints accounted for)
stephalf <- 0.215     # ex. 0.215 is around 200 full beam footprints, 0.215 seconds on each side of the "main" footprint. Change the value accordingly to the number of footprints you want. ex.: 0.215 seconds on each side of the "main" footprint

# Set the approach to be used
approach <- "singlebeam"  # approach singlebeam uses only neighboring footprints of the same beam 
approach <- "allbeams"    # approach allbeams uses the neighboring footprints of all beams of the dataset
approach <- "twobeams"    # approach twobeams uses the neighboring footprints of same laser unit beams of the dataset


#------------------------------------------
# Data preparation
#------------------------------------------

## Filter only good quality footprints
gedidata<- subset(gedidata,qlt_f=="01")           # quality_flag
gedidata<- subset(gedidata,dgrd_f=="00")          # degrage_flag
gedidata <- subset(gedidata,power == "full")      # full power beam : "BEAM0011",  "BEAM0010", "BEAM0001", "BEAM0000"

gedidata <- gedidata %>%
  mutate(lon = unlist(map(gedidata$geometry,1)),
         lat = unlist(map(gedidata$geometry,2)))

## Keep only essential data
gedidata <- gedidata %>% 
  dplyr::select(shotnumber,lon,lat,orbit,delta_time,nom_beam,elevG)    # elevG : elev_lowestmode corrected for height to fit coordinate system


gedidata <- as.data.frame(gedidata)
gedidata$geometry <- NULL

orbits <- unique(gedidata$orbit)


## Create search window
SearchSeq <- seq(-SearchDist, SearchDist, SearchStep)
SearchSeq.df <- merge(SearchSeq, SearchSeq, by=NULL)
colnames(SearchSeq.df) <- c("x_offset","y_offset")

nbextracted <- nrow(SearchSeq.df)


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
  
  chemin <- paste0(StatFold,sht,"rasterRStat2.tif")
  writeRaster(r,chemin,overwrite=T)
  
  output <- paste0(StatFold,sht,"accumulation2.tif")
  
  # apply flow accumulation FD8  
  whitebox::wbt_fd8_flow_accumulation(
    chemin,
    output=output,
    out_type = "cells",
    exponent = 1.1,
    threshold = NULL,
    log = FALSE,
    clip = FALSE,
    wd = StatFold,
    verbose_mode = FALSE
  )
  accum <- paste0(StatFold,sht,"accumulation2.tif")
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
    x_pond=SearchStep * round((weighted.mean(maxCell[,1],accum.raster[IdcellMax]))/SearchStep) # round to the same as the SearchStep
    y_pond=SearchStep * round((weighted.mean(maxCell[,2],accum.raster[IdcellMax]))/SearchStep)
    maxCell.df <- data.frame(x_pond,y_pond)
  }
  
  colnames(maxCell.df) <- c("x_offset","y_offset")
  accum.raster[IdcellMax]
  accum.raster <- setMinMax(accum.raster)
  IdcellMaxF <- maxValue(accum.raster)
  maxCell.df2 <- (cbind(maxCell.df,IdcellMaxF))
  
  return(maxCell.df2)
}  # end fonction



#---------------------------
# Parallel
#---------------------------
## extract elevation values orbit by orbit
# first we extract the elevation value at the initial position of each footprint. 
# Only if this value is in the DTM, so only if the footprint overlays the DTM raster, 
# then we also extract the elevation values of the surrounding search window

library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

z0 <- Sys.time()

processOrbit <- function(orb) {
  library(sp)
  library(dplyr)
  library(raster)
  library(terra)
  library(bit64)
  library(plyr)
  library(whitebox)
  library(ModelMetrics)

  MNTsmooth <- "\\RGEAlti\\RGEAlti_D11bufferMNTmean.tif"
  #with terra package
  MNTsmooth <- rast(MNTsmooth)
  crs(MNTsmooth)  <- "epsg:2154"
  
  #---------------------------
  # Extract elevation values
  #---------------------------
  
  allpoints <- data.frame()
  
  print(orb)
  todo <- gedidata %>% 
    filter(orbit == orb)
  SP_ij <- SpatialPoints(cbind(todo$lon, todo$lat))
  SP_ij <- vect(SP_ij) #with terra
  todo$elev00 <- terra::extract(MNTsmooth, SP_ij)[2] # extraction at initial position with terra
  #todo$elev <- raster::extract(MNTsmooth, SP_ij) # extraction at initial position with raster
  todo$elev00 <- unlist(todo$elev00)
  
  todo2 <- subset(todo,!is.na(elev00)) # if the footprints overlays the raster
  allpoints <- merge(todo2, SearchSeq.df, by=NULL) #then the search window is defined
  rm(todo)
  rm(todo2)
  gc()
  allpoints$latrecalc <- allpoints$lat + allpoints$y_offset 
  allpoints$lonrecalc  <- allpoints$lon + allpoints$x_offset
  SP_ij <- SpatialPoints(cbind(allpoints$lonrecalc, allpoints$latrecalc))
  SP_ij <- vect(SP_ij) #with terra
  allpoints$elev <- terra::extract(MNTsmooth, SP_ij)[2] # extraction at all positions defined in search window settings
  allpoints$elev <- unlist(allpoints$elev)
  
  allpoints2 <- allpoints %>%  # only keep footprint if a value could be extracted for all positions of the search window
    group_by(shotnumber) %>%
    filter(!any(is.na(elev))) %>%
    dplyr::summarise(count      = n())
  allpoints2 <- subset(allpoints2,count==nbextracted)     # nbextracted ==  2601 with SearchDist 50 and Searchstep 2
  allpoints <- allpoints[allpoints$shotnumber %in% allpoints2$shotnumber,]
  
  
  saveRDS(allpoints,file = paste0("/GEDI/VosgesD11b/",outputfilename,orb,".rds"))
  gc()
   
  # join gedi data and DEMref allpoints
  #gedidataAP <- left_join(allpoints, gedidata, by=c("shotnumber"))
  gedidataAP <- allpoints
  rm(allpoints)
  #get difference between DEMref elevation (elev) and GEDI elev_lowest (elevG)
  gedidataAP$diff <- gedidataAP$elev - gedidataAP$elevG
  
  #delete if diff > 100
  gedidataAP <- gedidataAP %>%
    filter(abs(diff) <=100)
  
  gedidataAP2 <- gedidataAP %>%  # delete if more than 100m difference
    group_by(shotnumber) %>%
    filter(!any(is.na(elev))) %>%
    dplyr::summarise(count      = n())
  gedidataAP2 <- subset(gedidataAP2,count==nbextracted) #nbextracted ==  2601 with 50 and 2
  gedidataAP <- gedidataAP[gedidataAP$shotnumber %in% gedidataAP2$shotnumber,]
  
  
  # General mean stat by orbit (with initial position (x_offset = 0 and y_offset = 0)
  Stats0 <- gedidataAP %>%
    filter (x_offset==0 & y_offset==0) %>%
    group_by(orbit) %>%
    dplyr::summarise(footprint_nb      = n(),
                     shift_x= mean(x_offset),
                     shift_y= mean(y_offset),
                     dtm_mean_error = mean(diff),
                     dtm_mean_abs_error = mean(abs(diff)), #diff = elev-elevG
                     dtm_mean_corr = cor(elev,elevG),
                     dtm_mean_rmse = ModelMetrics::rmse(elev,elevG)
    )
  
  # General mean stat by orbit and tested position to test by total orbit (without use of neigh_steptime)
  getOpt <- gedidataAP %>%
    group_by(orbit,x_offset,y_offset) %>%
    dplyr::summarise(footprint_nb      = n(),
                     shift_x= mean(x_offset),
                     shift_y= mean(y_offset),
                     Err = mean(diff),
                     AbsErr = mean(abs(diff)),
                     Corr = cor(elev,elevG),
                     RMSE = ModelMetrics::rmse(elev,elevG)
    )
  
  #------------------------------------------
  # Preparation of Flow accumulation execution 
  #------------------------------------------
  
  # get general optimal position for all footprints combined
  getOpt_accum <- getOpt %>%
    group_by(orbit) %>%
    do(flowaccum(.,StatFold=StatFold,criteria = "bary", var='AbsErr',sht=orb))
  
  optim_accum <- merge(gedidataAP, getOpt_accum [,c("orbit","x_offset","y_offset")], by=c("orbit","x_offset","y_offset"))
  
  # keep only needed variables
  gedidata_tile <- gedidataAP %>% 
    dplyr::select(orbit,nom_beam,shotnumber,delta_time,elevG,x_offset,y_offset,latrecalc,lonrecalc,elev,diff)
  
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
    if (approach == "twobeams"){
      if (nom_beamftp == "BEAM0101" | nom_beamftp == "BEAM0110"){
        gedidata_tilespec <- gedidata_tile %>% 
          filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          filter(nom_beam == "BEAM0101" | nom_beam =="BEAM0110")
      }
      if (nom_beamftp == "BEAM1000" | nom_beamftp == "BEAM1011"){
        gedidata_tilespec <- gedidata_tile %>% 
          filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          filter(nom_beam == "BEAM1000" | nom_beam =="BEAM1011")
      }
      
    }
    
    # calculate statistics for each position in search window
    getOpt_tilespec <- gedidata_tilespec %>% #for cluster
      group_by(orbit,x_offset,y_offset) %>%
      dplyr::summarise(footprint_nb      = n(),
                       shift_x= mean(x_offset),
                       shift_y= mean(y_offset),
                       Err = mean(diff),
                       AbsErr = mean(abs(diff)),
                       Corr = cor(elev,elevG),
                       RMSE = ModelMetrics::rmse(elev,elevG)
      )
    
    #if at least X footprints in the search window
    if (getOpt_tilespec$footprint_nb[1] >= 1){ # ex.: 1 is executed for all footprints, may be changed to higher value (ex. 30) to not execute code for footprints with few neighbours.
      
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
      
      # # save a backup every 1000 footprints
      # if (footprint%%1000 == 0){
      #   filetosave <- paste0(StatFold,"GeoGEDIfootprintshifts215BaryMae_",footprint,".rds")
      #   saveRDS(footprintshift0215barymae,filetosave)
      # }
    }
  }
  
  # save the final file
  saveRDS(footprintshift0215barymae,file = paste0("\\GEDI\\VosgesD11b\\GeoGEDIfootprintshiftsD11_",orb,".rds"))
          
          
  gedidataGeoGEDI <- merge(gedidataAP,footprintshift0215barymae[c("sht_ftp","x_offset","y_offset","IdcellMaxF","footprint_nb","Err","AbsErr","Corr","RMSE")],by.x=c("shotnumber","x_offset","y_offset"),by.y=c("sht_ftp","x_offset","y_offset"))
  saveRDS(gedidataGeoGEDI,paste0("\\GEDI\\VosgesD11b\\GeoGEDIfootprintD11_",orb,".rds"))
  return(gedidataGeoGEDI)
}

results <- foreach(orb = orbits, .packages = c("dplyr", "terra")) %dopar% {
  processOrbit(orb)
}

stopCluster(cl)

finalResult <- do.call(rbind,results)


saveRDS(finalResult,"\\\GEDI\\VosgesD11b\\GeoGEDIfootprintD11_30_1_part1.rds")
