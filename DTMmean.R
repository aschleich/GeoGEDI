#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(terra))

library(tools)
library(stringr)
library(terra)

terraOptions(progress=0)

########## Mean smooth MNT in a 25m window ##########
radius <- 12.5
method <- "sum"
args <- commandArgs(trailingOnly = TRUE)

smooth <- function(path) {
  # Read input data and NA value
  ras <- terra::rast(path)
  na <- grep("NoData Value", describe(path), value = TRUE)
  na <- strsplit(na, "=")[[1]][2] |> as.numeric()

  # Prepare kernel and apply smoothing
  focal_weight_matrix <- terra::focalMat(ras, radius, "circle")
  if (method == "max") {
    focal_weight_matrix[which(focal_weight_matrix != 0)] <- 1
  }
  smooth_raster <- terra::focal(ras, focal_weight_matrix, fun = method)
  terra::crs(smooth_raster) <- terra::crs(ras)
  rm(ras)
  gc()

  # Write output with GDAL create options
  out_tif <- paste(file_path_sans_ext(path), "_smooth.tif", sep = "")
  copt <- c("COMPRESS=ZSTD", "ZLEVEL=1", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  terra::writeRaster(smooth_raster, out_tif, overwrite = TRUE, gdal = copt, NAflag = na)
}

for (f in args) {
  smooth(f)
}
