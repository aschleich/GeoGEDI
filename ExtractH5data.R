#!/usr/bin/env Rscript
#---------------------------
# Convert HH5 GEDI data to RDS file + shift z values using geoid
#---------------------------
library(dplyr)
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

quality_filter <- FALSE
col_names <- c(
  "beam_name", "shotnumber", "delta_time", "lat", "lon", "elev_gedi",
  "selectedalgo", "degrade_flag", "surface_flag", "quality_flag", "sensit",
  "elv_bs_f", "elv_b0_err", "dem_gedi", "selected_mode", "selected_mode_flag",
  "ls_treecover", "ms_nveg", "ms_nveg_sd", "ms_treecover", "ms_treecover_sd",
  "rh0", "rh1", "rh2", "rh3", "rh4", "rh5", "rh6", "rh7", "rh8", "rh9", "rh10",
  "rh11", "rh12", "rh13", "rh14", "rh15", "rh16", "rh17", "rh18", "rh19", "rh20",
  "rh21", "rh22", "rh23", "rh24", "rh25", "rh26", "rh27", "rh28", "rh29", "rh30",
  "rh31", "rh32", "rh33", "rh34", "rh35", "rh36", "rh37", "rh38", "rh39", "rh40",
  "rh41", "rh42", "rh43", "rh44", "rh45", "rh46", "rh47", "rh48", "rh49", "rh50",
  "rh51", "rh52", "rh53", "rh54", "rh55", "rh56", "rh57", "rh58", "rh59", "rh60",
  "rh61", "rh62", "rh63", "rh64", "rh65", "rh66", "rh67", "rh68", "rh69", "rh70",
  "rh71", "rh72", "rh73", "rh74", "rh75", "rh76", "rh77", "rh78", "rh79", "rh80",
  "rh81", "rh82", "rh83", "rh84", "rh85", "rh86", "rh87", "rh88", "rh89", "rh90",
  "rh91", "rh92", "rh93", "rh94", "rh95", "rh96", "rh97", "rh98", "rh99", "rh100",
  "energy_total", "rx_maxamp"
)


time_start <- Sys.time()
final_df <- data.frame()

# Loop on all files (each file one orbit)
for (file in files) {
  print(paste0("Processing ", file))

  h5_obj <- rhdf5::H5Fopen(file)
  gedi_data <- data.frame()

  # Loop all the beams of one file (=all 8 beams from one orbit)
  for (j in 1:8) {
    stru2 <- h5dump(h5_obj, recursive = 2, all = FALSE, load = FALSE)
    beam_name <- names(stru2)[j]

    get_dataset <- function(variable) {
      name <- paste0("h5_obj&'/", beam_name, variable, "'")
      #print(name)
      path <- eval(parse(text = name))
      if (variable == "/shot_number") {
        return(rhdf5::H5Dread(path, bit64conversion = "bit64"))
      } else {
        return(rhdf5::H5Dread(path))
      }
    }

    lat_lowestmode <- get_dataset("/lat_lowestmode")
    lon_lowestmode <- get_dataset("/lon_lowestmode")
    elev_lowestmode <- get_dataset("/elev_lowestmode")
    degrade_flag <- get_dataset("/degrade_flag")
    dem_gedi <- get_dataset("/digital_elevation_model")
    surface_flag <- get_dataset("/surface_flag")
    quality_flag <- get_dataset("/quality_flag")
    shotnumber <- get_dataset("/shot_number")
    sensit <- get_dataset("/sensitivity")
    delta_time <- get_dataset("/delta_time")
    elv_bs_f <- get_dataset("/elevation_bias_flag")
    elv_b0_err <- get_dataset("/elevation_bin0_error")
    selected_mode <- get_dataset("/selected_mode")
    selected_mode_flag <- get_dataset("/selected_mode_flag")
    selectedalgo <- get_dataset("/selected_algorithm")
    ls_treecover <- get_dataset("/land_cover_data/landsat_treecover")
    ms_nveg <- get_dataset("/land_cover_data/modis_nonvegetated")
    ms_nveg_sd <- get_dataset("/land_cover_data/modis_nonvegetated_sd")
    ms_treecover <- get_dataset("/land_cover_data/modis_treecover")
    ms_treecover_sd <- get_dataset("/land_cover_data/modis_treecover_sd")
    # rh98 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[99,] #98%
    # rh100 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[101,] #100%
    # rh0 <- H5Dread(eval(parse(text=paste0("h5_obj&'/",beam_name,"/rh")[1,] #0%
    for (i in 0:100) {
      assign(paste0("rh", i), get_dataset("/rh")[i + 1, ])
    }
    energy_total <- get_dataset("/energy_total")
    rx_maxamp <- get_dataset("/rx_assess/rx_maxamp")

    gedi_data_beam <- data.frame(
      beam_name, shotnumber, delta_time, lat_lowestmode, lon_lowestmode, elev_lowestmode,
      selectedalgo, degrade_flag, surface_flag, quality_flag, sensit,
      elv_bs_f, elv_b0_err, dem_gedi, selected_mode, selected_mode_flag,
      ls_treecover, ms_nveg, ms_nveg_sd, ms_treecover, ms_treecover_sd,
      rh0, rh1, rh2, rh3, rh4, rh5, rh6, rh7, rh8, rh9, rh10,
      rh11, rh12, rh13, rh14, rh15, rh16, rh17, rh18, rh19, rh20,
      rh21, rh22, rh23, rh24, rh25, rh26, rh27, rh28, rh29, rh30,
      rh31, rh32, rh33, rh34, rh35, rh36, rh37, rh38, rh39, rh40,
      rh41, rh42, rh43, rh44, rh45, rh46, rh47, rh48, rh49, rh50,
      rh51, rh52, rh53, rh54, rh55, rh56, rh57, rh58, rh59, rh60,
      rh61, rh62, rh63, rh64, rh65, rh66, rh67, rh68, rh69, rh70,
      rh71, rh72, rh73, rh74, rh75, rh76, rh77, rh78, rh79, rh80,
      rh81, rh82, rh83, rh84, rh85, rh86, rh87, rh88, rh89, rh90,
      rh91, rh92, rh93, rh94, rh95, rh96, rh97, rh98, rh99, rh100,
      energy_total, rx_maxamp
    )
    names(gedi_data_beam) <- col_names

    gedi_data <- rbind(gedi_data, gedi_data_beam)

  } # end loop on beams

  # Filter elevation
  # gedi_data$diff_elevgedi <- gedi_data$elev_gedi - gedi_data$dem_gedi
  # gedi_data <- gedi_data %>% filter(diff_elevgedi < 300)

  # From data.frame to SpatVector
  gedi_data <- terra::vect(gedi_data, geom = c("lon", "lat"), crs = "epsg:4326")

  # Transform height to altitude using geoid undulation grid
  # h(IGN69) equals h(WGS84) minus h(RAF or RAC)
  gedi_data$elevgeoid_grid <- terra::extract(geoid_grid, gedi_data)
  gedi_data$elevRef <- gedi_data$elev_gedi - gedi_data$elevgeoid_grid

  # Reproject data to Lambert
  gedi_data_proj <- project(gedi_data, paste0("epsg:", epsg_code))

  # Extract columns of X and Y coordinates
  gedi_data$X_L93 <- terra::geom(gedi_data_proj)[, 3]
  gedi_data$Y_L93 <- terra::geom(gedi_data_proj)[, 4]
  gedi_data_proj$X_wgs <- terra::geom(gedi_data)[, 3]
  gedi_data_proj$Y_wgs <- terra::geom(gedi_data)[, 4]

  # Transform SpatialPointsDataFrame to data.frame
  gedi_data_proj.df <- as.data.frame(gedi_data_proj)

  # Keep only quality data
  if (quality_filter) {
    gedi_data_proj.df_nett <- subset(gedi_data_proj.df, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01" | surface_flag == "01" & degrade_flag == "01" & quality_flag == "01")
    gedi_data_proj.df_nett_strict <- subset(gedi_data_proj.df, surface_flag == "01" & degrade_flag == "00" & quality_flag == "01")
  }

  if (dim(gedi_data_proj.df)[1] == 0) {
    print("Couldn't find good quality footprints in our bbox")
  } else {
    print("Found good quality footprints in our bbox")
    final_df <- rbind(final_df, gedi_data_proj.df)
  } # end loop rbind data to final dataframe

  h5closeAll()
} # end loop for file

time_end <- Sys.time()
print(paste0("Elpased time : ", toString(time_end - time_start)))

final_df$orbit <- substr(final_df$shotnumber, 1, 4)

# Count number of different final orbits
nb_ftpt <- final_df %>% summarise(Unique_Elements = n_distinct(orbit))
print(paste0("Processed ", toString(nb_ftpt)))

# Add variable power to define full power and half power beam
half_beam_names <- c("BEAM0011", "BEAM0010", "BEAM0001", "BEAM0000")
final_df$power <- ifelse(final_df$beam_name %in% half_beam_names, "half", "full")

# Save entire dataframe, and one for full beams only
terra::saveRDS(final_df, file = "gedi_data_all.rds")
final_df_full <- final_df %>% filter(power == "full")
terra::saveRDS(final_df_full, file = "gedi_data_full.rds")
