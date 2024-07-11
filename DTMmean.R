#!/usr/bin/env Rscript
library(stringr)
# Terra env
suppressPackageStartupMessages(library(terra))
library(terra)
terraOptions(progress=0)

smooth <- function(path) {
  library(terra)
  ######## Mean smooth DTM in a 25m window ##########
  dtm <- terra::rast(path)
  
  # Prepare output
  nodata <- grep("NoData Value", terra::describe(path), value = TRUE)
  nodata <- strsplit(na, "=")[[1]][2] |> as.numeric()
  outfile <- paste(file_path_sans_ext(path), "_smooth.tif", sep = "")
  gdalopt <- c("COMPRESS=ZSTD", "ZLEVEL=1", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  
  # Create kernel
  radius <- 12.5
  kernel <- terra::focalMat(dtm, radius, "circle")
  method <- "sum" # or "max
  if (method == "max") {
    kernel[which(kernel != 0)] <- 1
  }
  
  # Use faster C++ implementation of "terra::focal"
  terra::focalCpp(dtm, kernel, method, fill_value=nodata, filename=outfile, wopt=list(gdalopt))
}

for (f in commandArgs(trailingOnly = TRUE)) {
  smooth(f)
}
