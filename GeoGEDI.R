#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(bit64, warn.conflicts = FALSE))
library(tidyr)
library(dplyr, warn.conflicts = FALSE)
use_arrow <- suppressPackageStartupMessages(require("arrow", warn.conflicts = FALSE))

options(digits = 12, error = traceback)

sep <- "/"
if (.Platform$OS.type == "windows") {
  sep <- "\\"
}

#------------------------------------------
# Input files
#------------------------------------------
arguments <- commandArgs(trailingOnly = TRUE)
gedidata_path <- arguments[1]
dem_smooth_path <- arguments[2]

target_crs <- terra::crs(terra::rast(dem_smooth_path))

#------------------------------------------
# Algorithm parameters
#------------------------------------------
n_cores <- 1

# Search window
search_dist <- 50
search_step <- 2

# Time step size in seconds (on each side of the "main" footprint)
# This fixes the time laps to consider neighboring footprints
# Change the value accordingly to the number of footprints you want
step_half <- 0.215
# 0.215 is around 200 full beam footprints, 0.215 seconds on each side of the "main" footprint

# Number of time steps forward to keep in sliding window of beam offsets
# Higher values will increase RAM usage but may speed up computation
steps_forward <- 4  # 4 = more or less 1GB of RAM per orbit / job in parallel

# Approach to be used while selecting footprints
# "singlebeam" uses only neighboring footprints of the same beam
# "allbeams" uses the neighboring footprints of all beams
# "twobeams" uses the neighboring footprints of same laser unit beams
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

nb_offsets <- nrow(search_df)

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
# Select footprints based on shot time and 'approach' setting
#------------------------------------------
filter_footprints <- function(gedidata_ap, time_ftp, beam_nameftp) {
  # Set time to select neighboring footprints for cluster
  time_ftpmin <- time_ftp - step_half
  time_ftpmax <- time_ftp + step_half
  df_neighbours <- gedidata_ap %>% dplyr::filter(delta_time > time_ftpmin, delta_time <= time_ftpmax)

  # Define cluster : select neighboring footprints depending on approach
  if (approach == "singlebeam") {
    return(dplyr::filter(df_neighbours, beam_name == beam_nameftp))
  } else if (approach == "allbeams") {
    return(df_neighbours)
  } else if (approach == "twobeams") {
    if (beam_nameftp == "BEAM0101" || beam_nameftp == "BEAM0110") {
      df_neighbours <- dplyr::filter(df_neighbours, beam_name == "BEAM0101" | beam_name == "BEAM0110")
    }
    if (beam_nameftp == "BEAM1000" || beam_nameftp == "BEAM1011") {
      df_neighbours <- dplyr::filter(df_neighbours, beam_name == "BEAM1000" | beam_name == "BEAM1011")
    }
    if (beam_nameftp == "BEAM0000" || beam_nameftp == "BEAM0001") {
      df_neighbours <- dplyr::filter(df_neighbours, beam_name == "BEAM0000" | beam_name == "BEAM0001")
    }
    if (beam_nameftp == "BEAM0010" || beam_nameftp == "BEAM0011") {
      df_neighbours <- dplyr::filter(df_neighbours, beam_name == "BEAM0010" | beam_name == "BEAM0011")
    }
    if (length(unique(df_neighbours$beam_name)) == 2) {
      return(df_neighbours)
    } else {
      # message("Only one beam name left after filtering, cannot use the 'twobeams' approach.")
      return(NULL)
    }
  } else {
    quit("Wrong value for parameter 'approach', possible values are 'singlebeam', 'twobeams', 'allbeams'.")
  }
}

get_offsets <- function(df_todo, dem_smooth) {
  df_todo <- dplyr::cross_join(df_todo, search_df)

  df_todo$x_shifted <- df_todo$x + df_todo$x_offset
  df_todo$y_shifted <- df_todo$y + df_todo$y_offset

  # Extraction at all positions defined in search window settings
  shifted_points <- terra::vect(dplyr::select(df_todo, x_shifted, y_shifted), geom = c("x_shifted", "y_shifted"), crs = target_crs)
  df_todo$elev <- unlist(terra::extract(dem_smooth, shifted_points)[2])

  # Get diff between DEMref elevation (elev) and GEDI elev_lowest (elev_ngf)
  # elev_ngf : elev_lowestmode corrected with geoid / vertical CRS shift
  df_todo$diff <- df_todo$elev - df_todo$elev_ngf
  # Delete if diff > 100
  df_todo <- df_todo %>% dplyr::filter(abs(diff) <= 100)

  # Only keep footprint if a DEM value could be extracted for all positions
  df_summary <- df_todo %>%
    dplyr::group_by(shot_number) %>%
    dplyr::filter(!any(is.na(elev))) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count == nb_offsets)

  df_todo <- df_todo %>%
    dplyr::filter(shot_number %in% df_summary$shot_number) %>%
    dplyr::select(orbit, shot_number, delta_time, elev_ngf, x_offset, y_offset, x_shifted, y_shifted, elev, diff)

  return(df_todo)
}

#------------------------------------------
# Flow accumulation algorithm
#------------------------------------------
flowaccum <- function(df, accum_dir, criterion, variable, shot) {
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
    plot_bounds <- c(-search_dist, search_dist)
    terra::plot(errormap, main = basename, asp = 1, xlim = plot_bounds, ylim = plot_bounds)
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

  accum <- terra::rast(accum_path)

  if (error_plots) {
    basename <- paste0(accum_dir, sep, "mnt_flow_", shot, "_", orb, "_", variable)
    grDevices::png(filename = paste0(basename, ".png"))
    terra::plot(accum, main = basename, asp = 1, xlim = c(-search_dist, search_dist))
    grDevices::dev.off()
  }

  # Select final "optimal" pixel out of flow accumulation map
  if (criterion == "min") {
    accum <- terra::setMinMax(accum)
    cell <- base::which.max(terra::values(accum))
    cell_xy <- terra::xyFromCell(accum, cell)
    max_cell_df <- data.frame(cell_xy)
  } else if (criterion == "bary") {
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

#------------------------------------------
# Main function to loop with, read and process one file per orbit
#------------------------------------------
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

  # Stop here if the footprints does not intersect with DTM
  gedidata_ap <- dplyr::filter(gedidata_ap, !is.na(elev00))
  if (dim(gedidata_ap)[1] == 0) {
    return(NULL)
  }

  # Preparation of Flow accumulation execution
  #------------------------------------------
  # Order by shot time
  gedidata_ap <- dplyr::arrange(gedidata_ap, delta_time)

  # Count number of footprints
  nb_ftp <- gedidata_ap %>%
    dplyr::summarise(Unique_Elements = dplyr::n_distinct(shot_number)) %>%
    as.data.frame()
  nb_ftp <- nb_ftp[, 1]

  # Flow accumulation algorithm applied to each footprint
  #------------------------------------------

  # Sliding window df of neighbours with offsets
  df_offsets <- NULL
  # Results df: optimal offset for each footprint
  df_results <- NULL

  for (ftp_idx in 1:nb_ftp) {
    shot_number <- gedidata_ap[ftp_idx, ]$shot_number
    time_ftp <- gedidata_ap[ftp_idx, ]$delta_time
    beam_nameftp <- gedidata_ap[ftp_idx, ]$beam_name
    win_time_max <- time_ftp + step_half * (steps_forward + 1)

    df_neighbours <- filter_footprints(gedidata_ap, time_ftp, beam_nameftp)
    if (is.null(df_neighbours) || nrow(df_neighbours) == 0) {
      # In case df is empty or single beam name left but approach is twobeams
      next
    }

    last_beam <- df_neighbours[df_neighbours$delta_time == max(df_neighbours$delta_time), ]
    if (is.null(df_offsets)) {
      df_todo <- gedidata_ap %>% dplyr::filter(delta_time > time_ftp - step_half, delta_time <= win_time_max)
      df_offsets <- get_offsets(df_todo, dem_smooth)
    } else if (!(last_beam[1]$shot_number %in% df_offsets$shot_number)) {
      df_offsets <- dplyr::filter(df_offsets, delta_time > time_ftp - step_half)
      df_todo <- gedidata_ap %>% dplyr::filter(delta_time > time_ftp - step_half, delta_time <= win_time_max)
      df_todo <- dplyr::filter(df_todo, !(shot_number %in% df_offsets$shot_number))
      df_offsets <- rbind(df_offsets, get_offsets(df_todo, dem_smooth))
    }

    df_current_offsets <- filter(df_offsets, shot_number %in% df_neighbours$shot_number)
    if (is.null(df_current_offsets) || nrow(df_current_offsets) == 0) {
      next
    }

    # Calculate statistics for each position in search window
    gedidata_tile <- df_current_offsets %>%
      dplyr::group_by(orbit, x_offset, y_offset) %>%
      dplyr::summarise(
        .groups = "keep",
        footprint_nb = n(),
        AbsErr = mean(abs(diff)),
        #Err = mean(diff),
        #Corr = cor(elev, elev_ngf),
        RMSE = ModelMetrics::rmse(elev, elev_ngf)
      )

    # If at least X footprints in the search window
    # Use 1 to execute on all footprints, may be changed to higher value (ex. 30)
    # to avoid running code for footprints with few neighbours.
    if (gedidata_tile$footprint_nb[1] >= 1) {
      # Run flow accumulation algorithm to find best shift
      # Define bary or min criterion to select final optimal position
      # Define variable ex.: AbsErr to create error map on
      df_accum_bary <- gedidata_tile %>%
        dplyr::group_by(orbit) %>%
        dplyr::do(flowaccum(., accum_dir = accum_dir, criterion = "bary", variable = "AbsErr", shot = shot_number))

      df_accum_bary$shot_number <- shot_number
      # Join
      df_accum_tilespec_bary <- dplyr::inner_join(df_accum_bary, df_current_offsets, by = c("orbit", "shot_number", "x_offset", "y_offset"))
      df_accum_tilespec_bary <- dplyr::inner_join(df_accum_tilespec_bary, gedidata_tile, by = c("orbit", "x_offset", "y_offset"))
      df_results <- rbind(df_results, df_accum_tilespec_bary)
    }
  }

  # Save the final file
  df_results <- dplyr::inner_join(gedidata_ap, df_results, by = c("orbit", "shot_number", "delta_time", "elev_ngf"))
  rm(gedidata_ap)
  df_results <- dplyr::inner_join(df_results, gedidata_full, by = c("orbit", "beam_name", "shot_number", "delta_time", "x", "y", "elev_ngf"))
  rm(gedidata_full)

  if (use_arrow) {
    arrow::write_parquet(make_arrow_table(df_results), paste0("O", orb, "_shifted.parquet"))
  } else {
    saveRDS(df_results, paste0("O", orb, "_shifted.rds"))
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

#------------------------------------------
# Execute script
#------------------------------------------
if (n_cores == 1 || length(gedi_files) < n_cores) {
  # Single core
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

#warnings()
