#!/usr/bin/env Rscript
#---------------------------
# Convert HH5 GEDI data to RDS file + shift z values using geoid
#---------------------------
library(dplyr, warn.conflicts = FALSE)
library(terra)
library(rhdf5)

# Script arguments
arguments <- commandArgs(trailingOnly = TRUE)

datadir <- arguments[1]
files <- paste0(datadir, "/", dir(datadir, pattern = "*.h5", recursive = TRUE))

# Download IGN geoid undulation grid, converted from txt to GeoTIFF
# See PROJ-data github repository: https://github.com/OSGeo/PROJ-data/tree/master/fr_ign
geoid_grid_file <- arguments[2]

if (startsWith(basename(geoid_grid_file), "fr_ign_RAF")) {
  epsg_code <- "5698"
} else if (startsWith(basename(geoid_grid_file), "fr_ign_RAC")) {
  epsg_code <- "5699"
} else {
  print("Could not guess CRS based on IGN geoid grid name.")
  quit(status = 1)
}
geoid_grid <- terra::rast(geoid_grid_file)

max_elevation_diff <- 300
quality_filter <- FALSE
out_vector <- TRUE
variables <- c(
  "/lat_lowestmode",
  "/lon_lowestmode",
  "/elev_lowestmode",
  "/degrade_flag",
  "/digital_elevation_model",
  "/surface_flag",
  "/quality_flag",
  "/shot_number",
  "/sensitivity",
  "/delta_time",
  "/elevation_bias_flag",
  "/elevation_bin0_error",
  "/selected_mode",
  "/selected_mode_flag",
  "/selected_algorithm",
  "/land_cover_data/landsat_treecover",
  "/land_cover_data/modis_nonvegetated",
  "/land_cover_data/modis_nonvegetated_sd",
  "/land_cover_data/modis_treecover",
  "/land_cover_data/modis_treecover_sd",
  "/energy_total",
  "/rx_assess/rx_maxamp"
)


time_start <- Sys.time()
final_df <- data.frame()

# Loop on all files (each file one orbit)
for (file in files) {
  print(paste0("Processing ", file))

  tryCatch(
    {
      h5_obj <- rhdf5::H5Fopen(file)
      gedi_data <- data.frame()

      # Loop all the beams of one file (=all 8 beams from one orbit)
      for (j in 1:8) {
        stru2 <- h5dump(h5_obj, recursive = 2, all = FALSE, load = FALSE)
        beam_name <- names(stru2)[j]

        get_dataset <- function(variable) {
          name <- paste0("h5_obj&'/", beam_name, variable, "'")
          path <- eval(parse(text = name))
          if (variable == "/shot_number") {
            return(rhdf5::H5Dread(path, bit64conversion = "bit64"))
          } else {
            return(rhdf5::H5Dread(path))
          }
        }

        data <- lapply(variables, get_dataset)
        col_names <- unlist(lapply(variables, function(varname) tail(unlist(strsplit(varname, "/")), 1)))
        # rh98 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[99,] #98%
        # rh100 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[101,] #100%
        # rh0 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[1,] #0%
        for (i in 0:100) {
          data <- append(data, list(get_dataset("/rh")[i + 1, ]))
          col_names <- append(col_names, paste0("rh", i))
        }
        data <- c(beam_name, data)
        col_names <- c("beam_name", col_names)
        beam_data <- as.data.frame(data, col.names = col_names)
        gedi_data <- rbind(gedi_data, beam_data)
      } # end loop on beams

      # Filter elevation
      # gedi_data$diff_elevgedi <- gedi_data$elev_gedi - gedi_data$dem_gedi
      # gedi_data <- gedi_data %>% filter(diff_elevgedi < max_elevation_diff)

      # From data.frame to SpatVector
      gedi_data <- terra::vect(gedi_data, geom = c("lon_lowestmode", "lat_lowestmode"), crs = "epsg:4326")
      # Transform height to altitude using geoid undulation grid
      # h(IGN69) equals h(WGS84) minus h(RAF or RAC)
      geoid_crs <- terra::crs(geoid_grid)
      suppressWarnings({
        gedi_data$elevgeoid_grid <- terra::extract(geoid_grid, project(gedi_data, geoid_crs))
      })
      gedi_data$elevRef <- gedi_data$elev_gedi - gedi_data$elevgeoid_grid

      # Reproject data to Lambert
      gedi_data_proj <- terra::project(gedi_data, paste0("epsg:", epsg_code))

      # Extract columns of X and Y coordinates
      gedi_data$X_L93 <- terra::geom(gedi_data_proj)[, 3]
      gedi_data$Y_L93 <- terra::geom(gedi_data_proj)[, 4]
      gedi_data_proj$X_wgs <- terra::geom(gedi_data)[, 3]
      gedi_data_proj$Y_wgs <- terra::geom(gedi_data)[, 4]

      # Transform SpatialPointsDataFrame to data.frame
      gedi_data_proj <- as.data.frame(gedi_data_proj)

      # Keep only quality data
      if (quality_filter) {
        # gedi_data_proj <- subset(gedi_data_proj, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01" | surface_flag == "01" & degrade_flag == "01" & quality_flag == "01")
        gedi_data_proj <- subset(gedi_data_proj, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01")
      }

      if (dim(gedi_data_proj)[1] == 0) {
        warning("Couldn't find good quality samples in our bbox")
      } else {
        print(paste("Found", toString(dim(gedi_data_proj)[1]), "samples"))
        final_df <- rbind(final_df, gedi_data_proj)
      } # end loop rbind data to final dataframe

      h5closeAll()
    },
    error = function(e) {
      print("Failed to process ", file)
      print(e)
    }
  )
} # end loop for file

time_end <- Sys.time()

final_df$orbit <- substr(final_df$shot_number, 1, 4)

# Count number of different final orbits
nb_orbit <- final_df %>% summarise(Unique_Elements = n_distinct(orbit))
elapsed_time <- round(time_end - time_start, 2)
print(paste("Processed", toString(nb_orbit), "orbits in", elapsed_time, "s"))

# Add variable power to define full power and half power beam
half_beam_names <- c("BEAM0011", "BEAM0010", "BEAM0001", "BEAM0000")
final_df$power <- ifelse(final_df$beam_name %in% half_beam_names, "half", "full")

# Save entire dataframe, and one for full beams only
terra::saveRDS(final_df, file = "gedi_data_all.rds")
final_df_full <- final_df %>% filter(power == "full")
terra::saveRDS(final_df_full, file = "gedi_data_full.rds")

if (out_vector) {
  gdf <- terra::vect(final_df, geom = c("X_wgs", "Y_wgs"), crs = "EPSG:4326")
  writeVector(gdf, "gedi_data_all.gpkg", overwrite = TRUE)
}
