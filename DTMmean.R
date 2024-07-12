#!/usr/bin/env Rscript
library(tools)
library(stringr)
# Terra env
suppressPackageStartupMessages(library(terra))
library(terra)
terraOptions(progress=0)

smooth <- function(path) {
  ######## Mean smooth DTM in a 25m window ##########
  outfile <- paste(file_path_sans_ext(path), "_smooth.tif", sep = "")
  if (!file.exists(outfile)) {
    dtm <- terra::rast(path)
    # Prepare output
    nodata <- -99999
    copt <- c("COMPRESS=ZSTD", "ZLEVEL=1", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
    # Create kernel
    radius <- 12.5
    kernel <- terra::focalMat(dtm, radius, "circle")
    method <- "sum" # for focal mean
    if (method == "max") {
      kernel[which(kernel != 0)] <- 1
    }
    # Apply focal function and write output
    a <- terra::focal(dtm, w=kernel, fun=method, filename=outfile, wopt=list(gdal=copt, NAflag=nodata))
  }
}

for (f in commandArgs(trailingOnly = TRUE)) {
  smooth(f)
}
