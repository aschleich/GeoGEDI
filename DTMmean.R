#!/usr/bin/env Rscript
library(tools)
library(stringr)
# Terra env
suppressPackageStartupMessages(library(terra))
library(terra)
terraOptions(progress=0)

########## Mean smooth DTM in a 25m window ##########
radius <- 12.5
method <- "sum"

smooth <- function(path) {
  # Read input data and NA value
  dtm <- terra::rast(path)
  nodata <- grep("NoData Value", describe(path), value = TRUE)
  nodata <- strsplit(na, "=")[[1]][2] |> as.numeric()
  
  # Prepare kernel
  kernel <- terra::focalMat(dtm, radius, "circle")
  if (method == "max") {
    kernel[which(kernel != 0)] <- 1
  }
  
  # Apply fun and directly write output
  outfile <- paste(file_path_sans_ext(path), "_smooth.tif", sep = "")
  copt <- c("COMPRESS=ZSTD", "ZLEVEL=1", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  
  # Use faster focal CPP implementation
  terra::focalCpp(dtm, kernel, method, fill_value=nodata, filename=outfile, wopt=list(copt))
}

args <- commandArgs(trailingOnly = TRUE)
for (f in args) {
  smooth(f)
}
