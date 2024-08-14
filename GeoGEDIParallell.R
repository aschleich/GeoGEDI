#!/usr/bin/env Rscript
library(plyr)
library(dplyr, warn.conflicts = FALSE)
library(stringr)
library(tidyverse)
library(sf)
library(sp)
suppressPackageStartupMessages(library(bit64, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(terra))
library(terra, warn.conflicts = FALSE)
library(whitebox)
library(ModelMetrics, warn.conflicts = FALSE)

# Force R to print more decimals
options(digits = 12)

if (.Platform$OS.type == "windows") {
  sep <- "\\"
} else {
  sep <- "/"
}

#------------------------------------------
# Input files
#------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)
gedidata_path <- normalizePath(arguments[1])
dem_smooth_path <- normalizePath(arguments[2])

target_crs <- crs(rast(dem_smooth_path))
quality_filter <- FALSE

#------------------------------------------
# Algorithm parameters
#------------------------------------------

# For parallelization
n_cores <- 1

# Set Search window
search_dist <- 50
search_step <- 2

# Set the time step size in seconds (on each side of the "main" footprint)
# This fixes the time laps to consider neighboring footprints when creating cluster
# Change the value accordingly to the number of footprints you want
step_half <- 0.215
# 0.215 is around 200 full beam footprints, 0.215 seconds on each side of the "main" footprint

# Set the approach to be used
# singlebeam uses only neighboring footprints of the same beam
approach <- "singlebeam"
# allbeams uses the neighboring footprints of all beams of the dataset
approach <- "allbeams"
# twobeams uses the neighboring footprints of same laser unit beams of the dataset
approach <- "twobeams"

#------------------------------------------
# Outputs
#------------------------------------------
error_plots <- TRUE

accum_dir <- "accum"
results_dir <- "results"

setwd(dirname(gedidata_path))
for (d in c(accum_dir, results_dir)) {
  if (!dir.exists(d)) {
    dir.create(d)
  }
}

## Keep only essential data
gedidata <- readRDS(gedidata_path)
if (quality_filter) {
  gedidata <- subset(gedidata, quality_flag == "01")
  gedidata <- subset(gedidata, degrage_flag == "00")
  gedidata <- subset(gedidata, power == "full")
}

gedidata <- gedidata %>%
  dplyr::select(shot_number, x, y, orbit, delta_time, beam_name, elev_ngf)
# elev_ngf : elev_lowestmode corrected with geoid / vertical CRS shift

orbits <- unique(gedidata$orbit)

## Create search window
search_seq <- seq(-search_dist, search_dist, search_step)
search_df <- merge(search_seq, search_seq, by = NULL)
colnames(search_df) <- c("x_offset", "y_offset")

nb_extracted <- base::nrow(search_df)

#------------------------------------------
# Flow accumulation algorithm function
#------------------------------------------

flowaccum <- function(df, accum_dir, criteria, var, shot) {
  df_1 <- df
  df_1$x_offset <- as.numeric(df_1$x_offset)
  df_1$y_offset <- as.numeric(df_1$y_offset)

  # create Error Map
  sp::coordinates(df_1) <- ~ x_offset + y_offset
  df_2 <- df_1[, c(var)]
  sp::gridded(df_2) <- TRUE
  r <- terra::rast(df_2, crs = target_crs)

  orb <- df_1$orbit[1]

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_", shot, "_", orb, "_", var)
    png(filename = paste0(basename, ".png"))
    plot(r, main = basename, asp = 1, xlim = c(-50, 50), lim = c(-50, 50))
    dev.off()
  }

  path <- paste0(accum_dir, sep, shot, "_error_map.tif")
  terra::writeRaster(r, path, overwrite = TRUE)

  accum_path <- paste0(accum_dir, sep, shot, "_accumulation.tif")
  Sys.sleep(3)

  # Apply flow accumulation FD8
  whitebox::wbt_fd8_flow_accumulation(
    path,
    output = accum_path,
    out_type = "cells",
    exponent = 1.1,
    threshold = NULL,
    wd = getwd(),
    log = FALSE,
    clip = FALSE,
    verbose_mode = TRUE
  )
  accum <- rast(accum_path)

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_flow_", shot, "_", orb, "_", var)
    png(filename = paste0(basename, ".png"))
    plot(accum, main = basename, asp = 1, xlim = c(-50, 50))
    dev.off()
  }

  # Select final "optimal" pixel out of flow accumulation map
  if (criteria == "min") {
    accum <- setMinMax(accum)
    id_cell_max <- which.max(accum)
    max_cell <- xyFromCell(accum, id_cell_max)
    max_cell_df <- data.frame(max_cell)
  } else if (criteria == "bary") {
    quant <- quantile(accum, probs = 0.99)
    id_cell_max <- which.max((accum) > quant)
    max_cell <- xyFromCell(accum, id_cell_max)
    # Weighted average flowaccumulation, rounded to the same as the search_step
    x_pond <- search_step * round((weighted.mean(max_cell[, 1], accum[id_cell_max])) / search_step)
    y_pond <- search_step * round((weighted.mean(max_cell[, 2], accum[id_cell_max])) / search_step)
    max_cell_df <- data.frame(x_pond, y_pond)
  }

  colnames(max_cell_df) <- c("x_offset", "y_offset")
  accum[id_cell_max]
  accum <- setMinMax(accum)
  id_cell_max_final <- minmax(accum)[2]
  max_cell_df <- (cbind(max_cell_df, id_cell_max_final))

  return(max_cell_df)
}

process_orbit <- function(orb) {
  ## Extract elevation values orbit by orbit
  # First we extract the elevation value at the initial position of each footprint.
  # Only if this value is in the DTM, so only if the footprint overlays the DTM raster,
  # then we also extract the elevation values of the surrounding search window
  dem_smooth <- rast(dem_smooth_path)

  #---------------------------
  # Extract elevation values
  #---------------------------
  allpoints <- gedidata %>% filter(orbit == orb)
  geopoints <- vect(allpoints, geom = c("x", "y"), crs = target_crs)
  # Extraction at initial position with terra
  allpoints$elev00 <- unlist(terra::extract(dem_smooth, geopoints)[2])
  allpoints <- subset(allpoints, !is.na(elev00)) # if the footprints overlays the raster
  if (dim(allpoints)[1] == 0) {
    return(NULL)
  }
  allpoints <- merge(allpoints, search_df, by = NULL) # then the search window is defined

  allpoints$x_shifted <- allpoints$x + allpoints$x_offset
  allpoints$y_shifted <- allpoints$y + allpoints$y_offset

  geopoints <- vect(allpoints, geom = c("x_shifted", "y_shifted"), crs = target_crs) # with terra
  # Extraction at all positions defined in search window settings
  allpoints$elev <- unlist(terra::extract(dem_smooth, geopoints)[2])

  # Only keep footprint if a value could be extracted for all positions of the search window
  allpoints2 <- allpoints %>%
    group_by(shot_number) %>%
    filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n())
  allpoints2 <- subset(allpoints2, count == nb_extracted)
  allpoints <- allpoints[allpoints$shot_number %in% allpoints2$shot_number, ]

  saveRDS(allpoints, file = paste0(results_dir, sep, orb, ".rds"))

  # Join gedi data and DEMref allpoints
  # gedidata_ap <- left_join(allpoints, gedidata, by=c("shot_number"))
  gedidata_ap <- allpoints
  rm(allpoints)
  gc()
  # Get difference between DEMref elevation (elev) and GEDI elev_lowest (elev_ngf)
  gedidata_ap$diff <- gedidata_ap$elev - gedidata_ap$elev_ngf

  # Delete if diff > 100
  gedidata_ap <- gedidata_ap %>% filter(abs(diff) <= 100)

  gedidata_ap2 <- gedidata_ap %>%
    group_by(shot_number) %>%
    filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n())
  gedidata_ap2 <- subset(gedidata_ap2, count == nb_extracted)
  gedidata_ap <- gedidata_ap[gedidata_ap$shot_number %in% gedidata_ap2$shot_number, ]

  # General mean stat by orbit (with initial position (x_offset = 0 and y_offset = 0)
  df_0 <- gedidata_ap %>%
    filter(x_offset == 0 & y_offset == 0) %>%
    group_by(orbit) %>%
    dplyr::summarise(
      footprint_nb = n(),
      shift_x = mean(x_offset),
      shift_y = mean(y_offset),
      dtm_mean_error = mean(diff),
      dtm_mean_abs_error = mean(abs(diff)), # diff = elev-elev_ngf
      dtm_mean_corr = cor(elev, elev_ngf),
      dtm_mean_rmse = ModelMetrics::rmse(elev, elev_ngf)
    )

  # General mean stat by orbit and tested position to test by total orbit (without use of neigh_steptime)
  df <- gedidata_ap %>%
    group_by(orbit, x_offset, y_offset) %>%
    dplyr::summarise(
      footprint_nb = n(),
      shift_x = mean(x_offset),
      shift_y = mean(y_offset),
      Err = mean(diff),
      AbsErr = mean(abs(diff)),
      Corr = cor(elev, elev_ngf),
      RMSE = ModelMetrics::rmse(elev, elev_ngf)
    )

  #------------------------------------------
  # Preparation of Flow accumulation execution
  #------------------------------------------

  # Get general optimal position for all footprints combined
  df_accum <- df %>%
    group_by(orbit) %>%
    do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = orb))

  optim_accum <- merge(gedidata_ap, df_accum[, c("orbit", "x_offset", "y_offset")], by = c("orbit", "x_offset", "y_offset"))

  # keep only needed variables
  gedidata_tile <- gedidata_ap %>%
    dplyr::select(orbit, beam_name, shot_number, delta_time, elev_ngf, x_offset, y_offset, y_shifted, x_shifted, elev, diff)

  # order the dataframes
  gedidata_tile <- gedidata_tile[order(gedidata_tile$delta_time), , drop = FALSE]
  optim_accum <- optim_accum[order(optim_accum$delta_time), , drop = FALSE]

  # count number of footprints
  gedidata_nbftp <- gedidata_tile %>%
    summarise(Unique_Elements = n_distinct(shot_number))
  gedidata_nbftp <- gedidata_nbftp[, 1]


  #------------------------------------------
  # Flow accumulation algorithm execution
  #------------------------------------------

  ftp_shift_bary <- data.frame()

  time <- Sys.time()

  for (footprint in 1:gedidata_nbftp) {
    shot_ftp <- optim_accum[footprint, ]$shot_number
    time_ftp <- optim_accum[footprint, ]$delta_time
    beam_nameftp <- optim_accum[footprint, ]$beam_name

    # Set time to select neighboring footprints for cluster
    time_ftpmin <- time_ftp - step_half
    time_ftpmax <- time_ftp + step_half

    # Define cluster : select neighboring footprints for calculation depending on approach
    if (approach == "singlebeam") {
      gedidata_tilespec <- gedidata_tile %>%
        filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
        filter(beam_name == beam_nameftp)
    }
    if (approach == "allbeams") {
      gedidata_tilespec <- gedidata_tile %>%
        filter(delta_time > time_ftpmin, delta_time <= time_ftpmax)
    }
    if (approach == "twobeams") {
      if (beam_nameftp == "BEAM0101" | beam_nameftp == "BEAM0110") {
        gedidata_tilespec <- gedidata_tile %>%
          filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          filter(beam_name == "BEAM0101" | beam_name == "BEAM0110")
      }
      if (beam_nameftp == "BEAM1000" | beam_nameftp == "BEAM1011") {
        gedidata_tilespec <- gedidata_tile %>%
          filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          filter(beam_name == "BEAM1000" | beam_name == "BEAM1011")
      }
    }

    # Calculate statistics for each position in search window
    df_tilespec <- gedidata_tilespec %>% # for cluster
      group_by(orbit, x_offset, y_offset) %>%
      dplyr::summarise(
        footprint_nb = n(),
        shift_x = mean(x_offset),
        shift_y = mean(y_offset),
        Err = mean(diff),
        AbsErr = mean(abs(diff)),
        Corr = cor(elev, elev_ngf),
        RMSE = ModelMetrics::rmse(elev, elev_ngf)
      )

    # If at least X footprints in the search window
    if (df_tilespec$footprint_nb[1] >= 1) { # ex.: 1 is executed for all footprints, may be changed to higher value (ex. 30) to not execute code for footprints with few neighbours.

      # Run flow accumulation algorithm to find best shift
      # Define bary or min criteria to select final optimal position
      # Define variable ex.: AbsErr to create error map on
      df_accum_bary <- df_tilespec %>%
        group_by(orbit) %>%
        do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = shot_ftp))

      # Join
      optim_accum_tilespec_bary <- merge(gedidata_tilespec, df_accum_bary[, c("orbit", "x_offset", "y_offset")], by = c("orbit", "x_offset", "y_offset"))

      df_accum_tilespec_bary <- join(df_accum_bary, df_tilespec, type = "inner")
      df_accum_tilespec_bary$shot_ftp <- shot_ftp
      df_accum_tilespec_bary$time_ftp <- time_ftp

      ftp_shift_bary <- rbind(ftp_shift_bary, df_accum_tilespec_bary)

      # Save a backup every 1000 footprints
      # if (footprint%%1000 == 0){
      #   saveRDS(ftp_shift_bary, paste0(accum_dir,"footprint_shift_bary_backup_", footprint, ".rds"))
      # }
    }
  }

  # Save the final file
  # saveRDS(ftp_shift_bary, file = paste0(results_dir, sep, "GeoGEDI_footprint_shift_full_", orb, ".rds"))

  col_names <- c("shot_ftp", "x_offset", "y_offset", "id_cell_max_final", "footprint_nb", "Err", "AbsErr", "Corr", "RMSE")
  by_x <- c("shot_number", "x_offset", "y_offset")
  by_y <- c("shot_ftp", "x_offset", "y_offset")
  GeoGEDI_data <- merge(gedidata_ap, ftp_shift_bary[col_names], by.x = by_x, by.y = by_y)
  # saveRDS(GeoGEDI_data, paste0(results_dir, sep, "GeoGEDI_footprints_", orb, ".rds"))
  return(GeoGEDI_data)
}

# Parallel execution
if (n_cores == 1) {
  # Single core - simply lapply
  final_result <- do.call(rbind, lapply(orbits, process_orbit))
} else {
  library(foreach, warn.conflicts = FALSE)
  library(doParallel)

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  libs <- c("plyr", "dplyr", "terra", "ModelMetrics")
  results <- foreach(orb = orbits, .packages = libs) %dopar% {
    process_orbit(orb)
  }
  stopCluster(cl)
}

saveRDS(final_result, paste0(results_dir, sep, "GeoGEDI_footprints.rds"))
