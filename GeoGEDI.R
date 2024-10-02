#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(bit64, warn.conflicts = FALSE))
library(tidyr)
library(dplyr, warn.conflicts = FALSE)
use_arrow <- suppressPackageStartupMessages(require("arrow", warn.conflicts = FALSE))

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
dem_smooth_path <- normalizePath(arguments[2])

target_crs <- terra::crs(terra::rast(dem_smooth_path))

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
# approach <- "singlebeam"
# allbeams uses the neighboring footprints of all beams of the dataset
# approach <- "allbeams"
# twobeams uses the neighboring footprints of same laser unit beams of the dataset
approach <- "twobeams"

# Filter already applied in ExtractH5data.R
quality_filter <- FALSE

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

## Create search window
search_seq <- seq(-search_dist, search_dist, search_step)
search_df <- merge(search_seq, search_seq, by = NULL)
colnames(search_df) <- c("x_offset", "y_offset")

nb_extracted <- nrow(search_df)

if (use_arrow) {
  # Convert f64 columns to f32 in order to save disk space
  make_arrow_table <- function(dataframe) {
    float64_cols <- sapply(dataframe, is.double)
    schema_list <- lapply(names(dataframe), function(col_name) {
      if (float64_cols[[col_name]] && !col_name %in% c("delta_time", "lat", "lon")) {
        arrow::field(col_name, arrow::float32())
      } else {
        arrow::field(col_name, arrow::infer_type(dataframe[[col_name]]))
      }
    })
    return(arrow::Table$create(dataframe, schema = arrow::schema(schema_list)))
  }
}

#------------------------------------------
# Flow accumulation algorithm function
#------------------------------------------
flowaccum <- function(df, accum_dir, criteria, variable, shot) {
  orb <- df$orbit[1]
  # Create Error Map
  sp::coordinates(df) <- ~ x_offset + y_offset
  df_err <- df[, c(variable)]
  sp::gridded(df_err) <- TRUE
  errormap <- terra::rast(df_err)
  terra::crs(errormap) <- target_crs

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_", shot, "_", orb, "_", variable)
    png(filename = paste0(basename, ".png"))
    terra::plot(errormap, main = basename, asp = 1, xlim = c(-50, 50), ylim = c(-50, 50))
    dev.off()
  }

  errormap_path <- paste0(accum_dir, sep, shot, "_errormap.tif")
  terra::writeRaster(errormap, errormap_path, overwrite = TRUE)

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
  accum <- terra::rast(accum_path)

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_flow_", shot, "_", orb, "_", variable)
    grDevices::png(filename = paste0(basename, ".png"))
    terra::plot(accum, main = basename, asp = 1, xlim = c(-50, 50))
    grDevices::dev.off()
  }

  # Select final "optimal" pixel out of flow accumulation map
  if (criteria == "min") {
    accum <- terra::setMinMax(accum)
    cell <- base::which.max(terra::values(accum))
    cell_xy <- terra::xyFromCell(accum, cell)
    max_cell_df <- data.frame(cell_xy)
  } else if (criteria == "bary") {
    quant <- terra::quantile(terra::values(accum), probs = 0.99)
    cell <- base::which.max(terra::values(accum > quant))
    cell_xy <- terra::xyFromCell(accum, cell)
    # Weighted average flowaccumulation, rounded to the same as the search_step
    x_pond <- search_step * round((stats::weighted.mean(cell_xy[, 1], accum[cell])) / search_step)
    y_pond <- search_step * round((stats::weighted.mean(cell_xy[, 2], accum[cell])) / search_step)
    max_cell_df <- data.frame(x_pond, y_pond)
  }

  accum <- terra::setMinMax(accum)
  max_accum <- terra::minmax(accum)[2]
  max_cell_df <- cbind(max_cell_df, max_accum)
  colnames(max_cell_df) <- c("x_offset", "y_offset", "max_accum")

  unlink(c(errormap_path, accum_path))
  return(max_cell_df)
}

process_footprint <- function(footprint_idx, gedidata_tile, optim_accum) {
  shot_number <- optim_accum[footprint_idx, ]$shot_number
  time_ftp <- optim_accum[footprint_idx, ]$delta_time
  beam_nameftp <- optim_accum[footprint_idx, ]$beam_name

  # Set time to select neighboring footprints for cluster
  time_ftpmin <- time_ftp - step_half
  time_ftpmax <- time_ftp + step_half
  gedidata_tile <- gedidata_tile %>%
    dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax)

  # Define cluster : select neighboring footprints depending on approach
  if (approach == "singlebeam") {
    gedidata_tile <- gedidata_tile %>%
      dplyr::filter(beam_name == beam_nameftp)
  } else if (approach == "twobeams") {
    if (beam_nameftp == "BEAM0101" || beam_nameftp == "BEAM0110") {
      gedidata_tile <- gedidata_tile %>%
        dplyr::filter(beam_name == "BEAM0101" | beam_name == "BEAM0110")
    }
    if (beam_nameftp == "BEAM1000" || beam_nameftp == "BEAM1011") {
      gedidata_tile <- gedidata_tile %>%
        dplyr::filter(beam_name == "BEAM1000" | beam_name == "BEAM1011")
    }
    if (beam_nameftp == "BEAM0000" || beam_nameftp == "BEAM0001") {
      gedidata_tile <- gedidata_tile %>%
        dplyr::filter(beam_name == "BEAM0000" | beam_name == "BEAM0001")
    }
    if (beam_nameftp == "BEAM0010" || beam_nameftp == "BEAM0011") {
      gedidata_tile <- gedidata_tile %>%
        dplyr::filter(beam_name == "BEAM0010" | beam_name == "BEAM0011")
    }
    if (length(unique(gedidata_tile$beam_name)) == 1) {
      # message("Only one beam name left after filtering, cannot use the 'twobeams' approach.")
      return(NULL)
    }
  } else if (approach != "allbeams") {
    stop("Wrong value for parameter 'approach', possible values are 'singlebeam', 'twobeams', 'allbeams'.")
  }

  # Calculate statistics for each position in search window
  gedidata_tile <- gedidata_tile %>% # for cluster
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
  # Use 1 to execute on all footprints, may be changed to higher value (ex. 30) to not execute code for footprints with few neighbours.
  if (gedidata_tile$footprint_nb[1] >= 1) {
    # Run flow accumulation algorithm to find best shift
    # Define bary or min criteria to select final optimal position
    # Define variable ex.: AbsErr to create error map on
    df_accum_bary <- gedidata_tile %>%
      dplyr::group_by(orbit) %>%
      dplyr::do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = shot_number))

    # Join
    df_accum_tilespec_bary <- dplyr::inner_join(df_accum_bary, gedidata_tile, by = c("orbit", "x_offset", "y_offset"))
    df_accum_tilespec_bary$shot_number <- shot_number
    df_accum_tilespec_bary$time_ftp <- time_ftp

    # Save a backup every 1000 footprints
    # if (footprint%%1000 == 0){
    #   saveRDS(ftp_shift_bary, paste0(accum_dir,"footprint_shift_bary_backup_", footprint, ".rds"))
    # }

    return(df_accum_tilespec_bary)
  }
}

process_orbit <- function(gedidata_path) {
  if (tools::file_ext(gedidata_path) == "parquet") {
    ext <- "parquet"
    gedidata_full <- arrow::read_parquet(gedidata_path)
  } else {
    ext <- "rds"
    gedidata_full <- tibble::as_tibble(readRDS(gedidata_path))
  }

  orb <- substr(basename(gedidata_path), 2, 6)
  if (file.exists(paste0("O", orb, "_shifted.", ext))) {
    message(paste("File", basename(gedidata_path), "already processed."))
    return(NULL)
  } else {
    message(paste("Processing file", basename(gedidata_path)))
  }

  if (quality_filter) {
    gedidata_full <- gedidata_full %>%
      dplyr::filter(surface_flag == "01", degrade_flag == "00", quality_flag == "01")
    gedidata_full <- gedidata_full %>%
      dplyr::select(surface_flag, degrade_flag, quality_flag)
    if (nrow(gedidata_full) == 0) {
      message("No remaining data after quality filter.")
      return(NULL)
    }
  }

  ## Keep only essential data
  gedidata_ap <- gedidata_full %>%
    dplyr::select(shot_number, x, y, orbit, delta_time, beam_name, elev_ngf)

  if (length(unique(gedidata_ap$orbit)) > 1) {
    stop("Invalid input file with multiple orbits.")
  }

  #---------------------------
  # Extract elevation values
  #---------------------------
  # First we extract the elevation value at the initial position of each footprint.
  # Only if this value is in the DTM, so only if the footprint overlays the DTM raster,
  # then we also extract the elevation values of the surrounding search window
  dem_smooth <- terra::rast(dem_smooth_path)

  # Extraction at initial position with terra
  gedidata_geo <- terra::vect(dplyr::select(gedidata_ap, x, y), geom = c("x", "y"), crs = target_crs)
  gedidata_ap$elev00 <- unlist(terra::extract(dem_smooth, gedidata_geo)[2])
  rm(gedidata_geo)

  gedidata_ap <- dplyr::filter(gedidata_ap, !is.na(elev00))

  # Stop here if the footprints does not intersect with DTM
  if (dim(gedidata_ap)[1] == 0) {
    return(NULL)
  }

  # Else the search window is defined
  gedidata_ap <- dplyr::cross_join(gedidata_ap, search_df)

  gedidata_ap$x_shifted <- gedidata_ap$x + gedidata_ap$x_offset
  gedidata_ap$y_shifted <- gedidata_ap$y + gedidata_ap$y_offset

  # Extraction at all positions defined in search window settings
  shifted_points <- terra::vect(dplyr::select(gedidata_ap, x_shifted, y_shifted), geom = c("x_shifted", "y_shifted"), crs = target_crs)
  gedidata_ap$elev <- unlist(terra::extract(dem_smooth, shifted_points)[2])
  rm(shifted_points, dem_smooth)

  # Get difference between DEMref elevation (elev) and GEDI elev_lowest (elev_ngf)
  # elev_ngf : elev_lowestmode corrected with geoid / vertical CRS shift
  gedidata_ap$diff <- gedidata_ap$elev - gedidata_ap$elev_ngf

  # Delete if diff > 100
  gedidata_ap <- gedidata_ap %>% dplyr::filter(abs(diff) <= 100)

  # Only keep footprint if a value could be extracted for all positions of the search window
  gedidata_summary <- gedidata_ap %>%
    dplyr::group_by(shot_number) %>%
    dplyr::filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count == nb_extracted)

  gedidata_ap <- gedidata_ap %>%
    dplyr::filter(shot_number %in% gedidata_summary$shot_number)

  if (dim(gedidata_ap)[1] == 0) {
    return(NULL)
  }

  rm(gedidata_summary)
  # saveRDS(gedidata_ap, file = paste0(results_dir, sep, orb, ".rds"))

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
  # Prepare a smaller dataset with only needed variables
  gedidata_tile <- gedidata_ap %>%
    dplyr::select(orbit, beam_name, shot_number, delta_time, elev_ngf, x_offset, y_offset, x_shifted, y_shifted, elev, diff)

  # Get general optimal position for all footprints combined
  df_accum <- df_accum %>%
    dplyr::group_by(orbit) %>%
    dplyr::do(flowaccum(., accum_dir = accum_dir, criteria = "bary", var = "AbsErr", shot = orb)) %>%
    dplyr::select(orbit, x_offset, y_offset)

  optim_accum <- merge(gedidata_ap, df_accum, by = c("orbit", "x_offset", "y_offset"))
  rm(df_accum)

  # Order the dataframes
  gedidata_tile <- gedidata_tile %>%
    dplyr::arrange(delta_time)
  optim_accum <- optim_accum %>%
    dplyr::arrange(delta_time)

  # Count number of footprints
  nb_ftp <- gedidata_tile %>%
    dplyr::summarise(Unique_Elements = dplyr::n_distinct(shot_number)) %>%
    as.data.frame()
  nb_ftp <- nb_ftp[, 1]

  #------------------------------------------
  # Flow accumulation algorithm applied to each footprint
  #------------------------------------------
  gedidata_shifted <- do.call(rbind, lapply(1:nb_ftp, process_footprint, gedidata_tile, optim_accum))
  rm(gedidata_tile, optim_accum)
  if (is.null(gedidata_shifted)) {
    return(NULL)
  }

  # Save the final file
  gedidata_shifted <- gedidata_shifted %>%
    dplyr::select(orbit, shot_number, x_offset, y_offset, max_accum, footprint_nb, Err, AbsErr, Corr, RMSE)
  # saveRDS(gedidata_shifted, file = paste0(results_dir, sep, "GeoGEDI_footprint_shift_full_", orb, ".rds"))
  geogedi_data <- dplyr::inner_join(gedidata_ap, gedidata_shifted, by = c("orbit", "shot_number", "x_offset", "y_offset"))
  rm(gedidata_ap, gedidata_shifted)
  geogedi_data <- dplyr::inner_join(geogedi_data, gedidata_full, by = c("orbit", "shot_number", "x", "y", "delta_time", "beam_name", "elev_ngf"))
  rm(gedidata_full)

  if (use_arrow) {
    arrow::write_parquet(make_arrow_table(geogedi_data), paste0("O", orb, "_shifted.parquet"))
  } else {
    saveRDS(geogedi_data, paste0("O", orb, "_shifted.rds"))
  }
}

if (use_arrow) {
  pattern <- "*.parquet"
} else {
  pattern <- "*.rds"
}

if (dir.exists(gedidata_path)) {
  gedi_files <- paste0(gedidata_path, sep, dir(gedidata_path, pattern = pattern))
} else if (file.exists(gedidata_path)) {
  gedi_files <- c(gedidata_path)
}

# Single core - simply lapply
if (n_cores == 1) {
  for (orb in sort(gedi_files)) {
    process_orbit(orb)
  }
  # Parallel execution
} else {
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)

  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  libs <- c("tidyr", "plyr", "dplyr", "terra", "ModelMetrics")
  # To print exported variables (duplicated in memory), use argument ".verbose = TRUE"
  foreach(orb = gedi_files, .packages = libs) %dopar% {
    process_orbit(orb)
  }
  stopCluster(cluster)
}
