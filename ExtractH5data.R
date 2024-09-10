#!/usr/bin/env Rscript
#################################################################
# Convert HDF5 GEDI data to RDS or parquet file + shift z values using geoid
#################################################################
suppressPackageStartupMessages(library(bit64, warn.conflicts = FALSE))
library(dplyr, warn.conflicts = FALSE)

use_arrow <- suppressPackageStartupMessages(require("arrow", warn.conflicts = FALSE))

# Script arguments
arguments <- commandArgs(trailingOnly = TRUE)
datadir <- normalizePath(arguments[1])
geoid_grid_file <- normalizePath(arguments[2])
roi_file <- normalizePath(arguments[3])

# First download IGN's grid of geoid undulation, converted from txt to GeoTIFF
# See PROJ-data github repository: https://github.com/OSGeo/PROJ-data/tree/master/fr_ign
if (startsWith(basename(geoid_grid_file), "fr_ign_RAF")) {
  epsg_code <- "5698"
} else if (startsWith(basename(geoid_grid_file), "fr_ign_RAC")) {
  epsg_code <- "5699"
} else {
  print("Could not guess CRS based on geoid filename.")
  epsg_code <- readline(prompt = "Enter your EPSG number:")
}
target_crs <- paste0("epsg:", epsg_code)

quality_filter <- TRUE
out_vector <- FALSE
variables <- c(
  "/shot_number",
  "/lat_lowestmode",
  "/lon_lowestmode",
  "/elev_lowestmode",
  "/digital_elevation_model",
  "/surface_flag",
  "/quality_flag",
  "/degrade_flag",
  "/sensitivity",
  "/delta_time",
  "/elevation_bias_flag",
  "/elevation_bin0_error",
  "/selected_mode",
  "/selected_mode_flag",
  "/selected_algorithm",
  "/energy_total",
  "/rx_assess/rx_maxamp",
  "/land_cover_data/landsat_treecover",
  "/land_cover_data/modis_nonvegetated",
  "/land_cover_data/modis_nonvegetated_sd",
  "/land_cover_data/modis_treecover",
  "/land_cover_data/modis_treecover_sd"
)

if (.Platform$OS.type == "unix") {
  sep <- "/"
} else {
  sep <- "\\"
}

# setwd(datadir)
if (!dir.exists("trash")) {
  dir.create("trash")
}

process_beam <- function(beam_num, h5_obj) {
  half_beam_names <- c("BEAM0011", "BEAM0010", "BEAM0001", "BEAM0000")
  # Read datasets headers only
  stru2 <- rhdf5::h5dump(h5_obj, recursive = 2, all = FALSE, load = FALSE)
  beam_name <- names(stru2)[beam_num]

  get_dataset <- function(ds_path) {
    name <- paste0("h5_obj&'/", beam_name, ds_path, "'")
    path <- eval(parse(text = name))
    if (ds_path == "/shot_number") {
      return(rhdf5::H5Dread(path, bit64conversion = "bit64"))
    } else {
      return(rhdf5::H5Dread(path))
    }
  }

  data <- lapply(variables, get_dataset)
  col_names <- unlist(lapply(variables, function(varname) tail(unlist(strsplit(varname, "/")), 1)))
  for (i in 0:100) {
    data <- append(data, list(get_dataset("/rh")[i + 1, ]))
    col_names <- append(col_names, paste0("rh", i))
  }
  data <- c(beam_name, data)
  col_names <- c("beam_name", col_names)
  beam_df <- as.data.frame(data, col.names = col_names)
  beam_df$power <- ifelse(beam_df$beam_name %in% half_beam_names, "half", "full")
  return(beam_df)
}

process_orbit <- function(h5_file) {
  # Support for zipped H5 files
  is_tmp <- FALSE
  if (tools::file_ext(h5_file) == "zip") {
    tmp <- tempdir()
    unzip(h5_file, exdir = tmp)
    h5_file <- paste0(tmp, sep, tools::file_path_sans_ext(basename(h5_file)))
    is_tmp <- TRUE
  }

  print(paste0("Processing ", h5_file))
  gedi_df <- NULL

  tryCatch(
    {
      h5_obj <- rhdf5::H5Fopen(h5_file)
      gedi_df <- do.call(rbind, lapply(1:8, process_beam, h5_obj))
      rhdf5::h5closeAll()
    },
    # When H5 file is corrupted
    error = function(e) {
      if (grepl("HDF5. Object header. Can't open object.", e)) {
        # TODO: check error message
        print(paste0("Invalid H5 file: ", h5_file))
        if (!is_tmp) {
          print("Moving to 'trash' folder")
          file.copy(h5_file, paste0("trash", sep, basename(h5_file)))
          file.remove(h5_file)
        }
      } else {
        print(e)
      }
    }
  )
  if (is.null(gedi_df)) return(NULL)

  gedi_df$orbit <- substr(gedi_df$shot_number, 1, 5)
  geom_cols <- c("lon_lowestmode", "lat_lowestmode")

  # Convert dataframe to SpatVector
  gedi_df <- terra::vect(gedi_df, geom = geom_cols, crs = "epsg:4326")
  # Filter points in polygons if an ROI file was provided
  if (!is.na(roi_file)) {
    gedi_df <- terra::mask(gedi_df, terra::vect(roi_file))
  }

  # Transform height to altitude using geoid undulation grid
  geoid_crs <- terra::crs(geoid_grid)
  suppressWarnings({
    gedi_df$geoid_height <- terra::extract(geoid_grid, terra::project(gedi_df, geoid_crs))[2]
  })
  # h(IGN69) = h(WGS84) - h(RAF18)
  gedi_df$elev_ngf <- gedi_df$elev_lowestmode - gedi_df$geoid_height

  # Save lon/lat and x/y coords
  gedi_df$lon <- terra::geom(gedi_df)[, 3]
  gedi_df$lat <- terra::geom(gedi_df)[, 4]
  gedi_df <- terra::project(gedi_df, target_crs)
  gedi_df$x <- terra::geom(gedi_df)[, 3]
  gedi_df$y <- terra::geom(gedi_df)[, 4]

  # Back to data.frame
  gedi_df <- terra::values(gedi_df)

  # Keep only quality data
  n_samples <- dim(gedi_df)[1]
  if (quality_filter) {
    gedi_df <- terra::subset(gedi_df, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01")
  }
  # Clean tmp data in case of unzip
  if (is_tmp) {
    file.remove(h5_file)
  }

  orbit <- unique(gedi_df$orbit)[1]
  if (dim(gedi_df)[1] == 0) {
    print("Couldn't find any good quality samples within ROI")
  } else {
    print(paste0("Extracted ", dim(gedi_df)[1], "/", n_samples, " samples"))

    if (use_arrow) {
      arrow::write_parquet(gedi_df, paste0(orbit,".parquet"))
    } else {
      saveRDS(gedi_df, file = paste0("O", orbit, ".rds"))
    }
    if (out_vector) {
      gedi_df <- terra::vect(gedi_df, geom = c("x", "y"), crs = target_crs)
      terra::writeVector(gedi_df, paste0("O", orbit, ".gpkg"), overwrite = TRUE)
    }
  }
}

geoid_grid <- terra::rast(geoid_grid_file)
files <- paste0(datadir, sep, dir(datadir, pattern = "*.h5", recursive = TRUE))

for (orb in files) { process_orbit(orb) }
