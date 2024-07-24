#!/usr/bin/env Rscript
#################################################################
# Convert HDF5 GEDI data to RDS file + shift z values using geoid
#################################################################

suppressPackageStartupMessages(library(bit64, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(terra))
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(rhdf5)

quality_filter <- TRUE
out_vector <- TRUE
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


process_beam <- function(h5_obj, beam_num) {
  half_beam_names <- c("BEAM0011", "BEAM0010", "BEAM0001", "BEAM0000")
  # Read datasets headers only
  stru2 <- h5dump(h5_obj, recursive = 2, all = FALSE, load = FALSE)
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
  h5_obj <- rhdf5::H5Fopen(h5_file)
  fun <- partial(process_beam, h5_obj)
  gedi_df <- as.data.frame(do.call(rbind, lapply(1:8, fun)))
  h5closeAll()

  gedi_df$orbit <- substr(gedi_df$shot_number, 1, 4)
  geom_cols <- c("lon_lowestmode", "lat_lowestmode")
  gedi_df <- vect(gedi_df, geom = geom_cols, crs = "epsg:4326")

  # Transform height to altitude using geoid undulation grid
  geoid_crs <- crs(geoid_grid)
  suppressWarnings({
    gedi_df$geoid_height <- extract(geoid_grid, project(gedi_df, geoid_crs))[2]
  })
  # h(IGN69) = h(WGS84) - h(RAF18)
  gedi_df$elev_ngf <- gedi_df$elev_lowestmode - gedi_df$geoid_height

  # Save lon/lat and x/y coords
  gedi_df$lon <- geom(gedi_df)[, 3]
  gedi_df$lat <- geom(gedi_df)[, 4]
  gedi_df <- project(gedi_df, paste0("epsg:", epsg_code))
  gedi_df$x <- geom(gedi_df)[, 3]
  gedi_df$y <- geom(gedi_df)[, 4]

  gedi_df <- as.data.frame(gedi_df)

  # Keep only quality data
  n_samples <- dim(gedi_df)[1]
  if (quality_filter) {
    gedi_df <- subset(gedi_df, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01")
  }

  if (dim(gedi_df)[1] == 0) {
    print("Couldn't find any good quality samples")
  } else {
    print(paste0("Extracted ", dim(gedi_df)[1], "/", n_samples, " samples"))
  } # end loop rbind data to final dataframe

  return(gedi_df)
}

# Script arguments
arguments <- commandArgs(trailingOnly = TRUE)
datadir <- arguments[1]

files <- paste0(datadir, "/", dir(datadir, pattern = "*.h5", recursive = TRUE))
if (!dir.exists("trash")) {
  dir.create("trash")
}

# First download IGN's grid of geoid undulation, converted from txt to GeoTIFF
# See PROJ-data github repository: https://github.com/OSGeo/PROJ-data/tree/master/fr_ign
geoid_grid_file <- arguments[2]
if (startsWith(basename(geoid_grid_file), "fr_ign_RAF")) {
  epsg_code <- "5698"
} else if (startsWith(basename(geoid_grid_file), "fr_ign_RAC")) {
  epsg_code <- "5699"
} else {
  print("Could not guess CRS based on geoid filename.")
  epsg_code <- readline(prompt = "Enter your EPSG number:")
}
geoid_grid <- rast(geoid_grid_file)

# Main loop on all files (one h5 file per orbit)
final_df <- data.frame()

for (f in files) {
  print(paste0("Processing ", f))
  tryCatch(
    {
      final_df <- rbind(final_df, process_orbit(f))
    },
    # When H5 file is corrupted
    error = function(e) {
      print(e)
      if (grepl("HDF5. Object header. Can't open object.", e)) {
        # TODO: check error message
        print(paste0("Failed to process ", f, ", moving to trash folder"))
        file.copy(f, paste0("trash/", basename(f)))
        file.remove(f)
      }
    }
  )
}

n_orbit <- final_df %>% summarise(Unique_Elements = n_distinct(orbit))
print(paste("Exporting", n_orbit, "orbits..."))

saveRDS(final_df, file = "gedi_data.rds")
if (out_vector) {
  gdf <- vect(final_df, geom = c("x", "y"), crs = paste0("epsg:", epsg_code))
  writeVector(gdf, "gedi_data.gpkg", overwrite = TRUE)
}
