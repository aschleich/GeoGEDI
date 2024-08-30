#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = FALSE)

options(digits = 12, error = traceback)

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
gedidata_path <- "/home/vidlb/Projets/umr_tetis/geogedi/test/GEDI_data.rds"
dem_smooth_path <- normalizePath(arguments[2])
dem_smooth_path <- "/home/vidlb/Projets/umr_tetis/geogedi/MNT/corse.vrt"

target_crs <- terra::crs(terra::rast(dem_smooth_path))
quality_filter <- FALSE

#------------------------------------------
# Algorithm parameters
#------------------------------------------
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
# approach <- "allbeams"
# twobeams uses the neighboring footprints of same laser unit beams of the dataset
# approach <- "twobeams"

#------------------------------------------
# Outputs
#------------------------------------------
error_plots <- FALSE

accum_dir <- "accum"
results_dir <- "results"

# setwd(dirname(gedidata_path))
for (d in c(accum_dir, results_dir)) {
  if (!dir.exists(d)) {
    dir.create(d)
  }
}

gedidata <- terra::readRDS(gedidata_path)
## Keep only essential data
gedidata <- gedidata %>%
  dplyr::select(shot_number, x, y, orbit, delta_time, beam_name, elev_ngf)
# elev_ngf : elev_lowestmode corrected with geoid / vertical CRS shift

if (quality_filter) {
  gedidata <- terra::subset(gedidata, quality_flag == "01" & degrage_flag == "00" & power == "full")
}

## Create search window
search_seq <- seq(-search_dist, search_dist, search_step)
search_df <- terra::merge(search_seq, search_seq, by = NULL)
colnames(search_df) <- c("x_offset", "y_offset")

nb_extracted <- base::nrow(search_df)

#------------------------------------------
# Flow accumulation algorithm function
#------------------------------------------
flowaccum <- function(df, accum_dir, criteria, var, shot) {
  # Create Error Map
  sp::coordinates(df) <- ~ x_offset + y_offset
  df_err <- df[, c(var)]
  sp::gridded(df_err) <- TRUE
  error_map <- terra::rast(df_err)
  terra::crs(error_map) <- target_crs

  orb <- df$orbit[1]

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_", shot, "_", orb, "_", var)
    png(filename = paste0(basename, ".png"))
    terra::plot(error_map, main = basename, asp = 1, xlim = c(-50, 50), ylim = c(-50, 50))
    dev.off()
  }

  errormap_path <- paste0(accum_dir, sep, shot, "_error_map.tif")
  terra::writeRaster(error_map, errormap_path, overwrite = TRUE)

  accum_path <- paste0(accum_dir, sep, shot, "_accumulation.tif")

  # Apply flow accumulation FD8
  whitebox::wbt_fd8_flow_accumulation(
    errormap_path,
    output = accum_path,
    out_type = "cells",
    exponent = 1.1,
    threshold = NULL,
    wd = getwd(),
    log = FALSE,
    clip = FALSE,
    compress_rasters = TRUE,
    verbose_mode = TRUE
  )
  # ! Use raster instead of terra because the quantile function doesn't work !
  accum <- raster::raster(accum_path)

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_flow_", shot, "_", orb, "_", var)
    grDevices::png(filename = paste0(basename, ".png"))
    terra::plot(accum, main = basename, asp = 1, xlim = c(-50, 50))
    grDevices::dev.off()
  }

  # Select final "optimal" pixel out of flow accumulation map
  if (criteria == "min") {
    accum <- raster::setMinMax(accum)
    id_cell_max <- raster::which.max(accum)
    max_cell <- raster::xyFromCell(accum, id_cell_max)
    max_cell_df <- data.frame(max_cell)
  } else if (criteria == "bary") {
    quant <- raster::quantile(accum, probs = 0.99)
    id_cell_max <- raster::which.max((accum) > quant)
    max_cell <- raster::xyFromCell(accum, id_cell_max)
    # Weighted average flowaccumulation, rounded to the same as the search_step
    x_pond <- search_step * round((weighted.mean(max_cell[, 1], accum[id_cell_max])) / search_step)
    y_pond <- search_step * round((weighted.mean(max_cell[, 2], accum[id_cell_max])) / search_step)
    max_cell_df <- data.frame(x_pond, y_pond)
  }

  colnames(max_cell_df) <- c("x_offset", "y_offset")
  accum[id_cell_max]
  accum <- raster::setMinMax(accum)
  id_cell_max_final <- raster::maxValue(accum)
  max_cell_df <- (cbind(max_cell_df, id_cell_max_final))

  return(max_cell_df)
}

process_footprint <- function(footprint_idx, gedidata_tile, optim_accum) {
    shot_ftp <- optim_accum[footprint_idx, ]$shot_number
    time_ftp <- optim_accum[footprint_idx, ]$delta_time
    beam_nameftp <- optim_accum[footprint_idx, ]$beam_name

    # Set time to select neighboring footprints for cluster
    time_ftpmin <- time_ftp - step_half
    time_ftpmax <- time_ftp + step_half

    # Define cluster : select neighboring footprints for calculation depending on approach
    if (approach == "singlebeam") {
      gedidata_tilespec <- gedidata_tile %>%
        dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
        dplyr::filter(beam_name == beam_nameftp)
    }
    if (approach == "allbeams") {
      gedidata_tilespec <- gedidata_tile %>%
        dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax)
    }
    if (approach == "twobeams") {
      if (beam_nameftp == "BEAM0101" || beam_nameftp == "BEAM0110") {
        gedidata_tilespec <- gedidata_tile %>%
          dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          dplyr::filter(beam_name == "BEAM0101" | beam_name == "BEAM0110")
      }
      if (beam_nameftp == "BEAM1000" || beam_nameftp == "BEAM1011") {
        gedidata_tilespec <- gedidata_tile %>%
          dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax) %>%
          dplyr::filter(beam_name == "BEAM1000" | beam_name == "BEAM1011")
      }
    }

    # Calculate statistics for each position in search window
    df_tilespec <- gedidata_tilespec %>% # for cluster
      dplyr::group_by(orbit, x_offset, y_offset) %>%
      dplyr::summarise(
        .groups = "keep",
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
        dplyr::group_by(orbit) %>%
        dplyr::do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = shot_ftp))

      # Join
      df_accum_tilespec_bary <- plyr::join(df_accum_bary, df_tilespec, type = "inner", by = c("orbit", "x_offset", "y_offset"))
      df_accum_tilespec_bary$shot_ftp <- shot_ftp
      df_accum_tilespec_bary$time_ftp <- time_ftp

      # Save a backup every 1000 footprints
      # if (footprint%%1000 == 0){
      #   saveRDS(ftp_shift_bary, paste0(accum_dir,"footprint_shift_bary_backup_", footprint, ".rds"))
      # }

      return(df_accum_tilespec_bary)
    }
}

process_orbit <- function(gedidata_ap) {
  ## Extract elevation values orbit by orbit
  orb <- unique(gedidata_ap$orbit)
  #print(paste("Processing orbit", orb))

  #---------------------------
  # Extract elevation values
  #---------------------------
  # First we extract the elevation value at the initial position of each footprint.
  # Only if this value is in the DTM, so only if the footprint overlays the DTM raster,
  # then we also extract the elevation values of the surrounding search window
  dem_smooth <- terra::rast(dem_smooth_path)

  # Extraction at initial position with terra
  gedidata_geo <- terra::vect(gedidata_ap, geom = c("x", "y"), crs = target_crs)
  gedidata_ap$elev00 <- unlist(terra::extract(dem_smooth, gedidata_geo)[2])
  rm(gedidata_geo)
  # drop = TRUE -> back to data.frame
  gedidata_ap <- terra::subset(gedidata_ap, !is.na(gedidata_ap$elev00), drop = TRUE)
  # Stop here if the footprints does not intersect with DTM
  if (dim(gedidata_ap)[1] == 0) {
    return(NULL)
  }
  # Else the search window is defined
  gedidata_ap <- base::merge(gedidata_ap, search_df, by = NULL)

  gedidata_ap$x_shifted <- gedidata_ap$x + gedidata_ap$x_offset
  gedidata_ap$y_shifted <- gedidata_ap$y + gedidata_ap$y_offset

  # Extraction at all positions defined in search window settings
  shifted_points <- terra::vect(gedidata_ap, geom = c("x_shifted", "y_shifted"), crs = target_crs)
  gedidata_ap$elev <- unlist(terra::extract(dem_smooth, shifted_points)[2])
  rm(shifted_points)

  # Only keep footprint if a value could be extracted for all positions of the search window
  gedidata_summary <- gedidata_ap %>%
    dplyr::group_by(shot_number) %>%
    dplyr::filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n())
  gedidata_summary <- base::subset(gedidata_summary, count == nb_extracted)
  gedidata_ap <- gedidata_ap[gedidata_ap$shot_number %in% gedidata_summary$shot_number, ]

  # saveRDS(gedidata_ap, file = paste0(results_dir, sep, orb, ".rds"))

  # Join gedi data and DEMref gedidata_ap
  # gedidata_ap <- dplyr::left_join(gedidata_ap, gedidata, by=c("shot_number"))
  # Get difference between DEMref elevation (elev) and GEDI elev_lowest (elev_ngf)
  gedidata_ap$diff <- gedidata_ap$elev - gedidata_ap$elev_ngf

  # Delete if diff > 100
  gedidata_ap <- gedidata_ap %>% dplyr::filter(abs(diff) <= 100)

  gedidata_summary <- gedidata_ap %>%
    dplyr::group_by(shot_number) %>%
    dplyr::filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n())
  gedidata_summary <- base::subset(gedidata_summary, count == nb_extracted)
  gedidata_ap <- gedidata_ap[gedidata_ap$shot_number %in% gedidata_summary$shot_number, ]
  rm(gedidata_summary)

  # General mean stat by orbit (with initial position (x_offset = 0 and y_offset = 0)
  # df_0 <- gedidata_ap %>%
  #   dplyr::filter(x_offset == 0 & y_offset == 0) %>%
  #   dplyr::group_by(orbit) %>%
  #   dplyr::summarise(
  #     .groups = "keep",
  #     footprint_nb = n(),
  #     shift_x = mean(x_offset),
  #     shift_y = mean(y_offset),
  #     dtm_mean_error = mean(diff),
  #     dtm_mean_abs_error = mean(abs(diff)), # diff = elev-elev_ngf
  #     dtm_mean_corr = cor(elev, elev_ngf),
  #     dtm_mean_rmse = ModelMetrics::rmse(elev, elev_ngf)
  #   )

  # General mean stat by orbit and tested position to test by total orbit (without use of neigh_steptime)
  df_accum <- gedidata_ap %>%
    dplyr::group_by(orbit, x_offset, y_offset) %>%
    dplyr::summarise(
      .groups = "keep",
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
  df_accum <- df_accum %>%
    dplyr::group_by(orbit) %>%
    dplyr::do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = orb))

  optim_accum <- base::merge(gedidata_ap, df_accum[, c("orbit", "x_offset", "y_offset")], by = c("orbit", "x_offset", "y_offset"))

  # Keep only needed variables
  gedidata_tile <- gedidata_ap %>%
    dplyr::select(orbit, beam_name, shot_number, delta_time, elev_ngf, x_offset, y_offset, x_shifted, y_shifted, elev, diff)

  # Order the dataframes
  gedidata_tile <- gedidata_tile[base::order(gedidata_tile$delta_time), , drop = FALSE]
  optim_accum <- optim_accum[base::order(optim_accum$delta_time), , drop = FALSE]

  # Count number of footprints
  nb_ftp <- gedidata_tile %>%
    dplyr::summarise(Unique_Elements = dplyr::n_distinct(shot_number))
  nb_ftp <- nb_ftp[, 1]

  #------------------------------------------
  # Flow accumulation algorithm applied to each footprint
  #------------------------------------------
  ftp_shift_bary <- do.call(rbind, lapply(1:nb_ftp, process_footprint, gedidata_tile, optim_accum))

  # Save the final file
  # saveRDS(ftp_shift_bary, file = paste0(results_dir, sep, "GeoGEDI_footprint_shift_full_", orb, ".rds"))

  col_names <- c("shot_ftp", "x_offset", "y_offset", "id_cell_max_final", "footprint_nb", "Err", "AbsErr", "Corr", "RMSE")
  by_x <- c("shot_number", "x_offset", "y_offset")
  by_y <- c("shot_ftp", "x_offset", "y_offset")
  geogedi_data <- terra::merge(gedidata_ap, ftp_shift_bary[col_names], by.x = by_x, by.y = by_y)

  return(geogedi_data)
}

# Group by orbit
gedidata <- gedidata %>%
  dplyr::group_by(orbit) %>%
  dplyr::group_split()

# Single core - simply lapply
if (n_cores == 1) {
  results <- do.call(rbind, lapply(gedidata, process_orbit))
# Parallel execution
} else {
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)

  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  libs <- c("plyr", "dplyr", "terra", "ModelMetrics")
  # To print exported variables (duplicated in memory), use argument ".verbose = TRUE"
  results <- foreach(group = gedidata, .packages = libs, .combine = rbind) %dopar% {
    process_orbit(group)
  }
  stopCluster(cluster)
}

saveRDS(results, "GeoGEDI_footprints.rds")
