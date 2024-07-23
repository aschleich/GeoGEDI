#!/usr/bin/env Rscript
library(tools)
library(stringr)
suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0)

smooth <- function(inpath, outpath) {
  ######## Mean smooth DTM in a 25m window ##########
  dtm <- terra::rast(inpath)
  # Prepare output
  nodata <- -99999
  # Create kernel
  radius <- 12.5
  kernel <- terra::focalMat(dtm, radius, "circle")
  method <- "sum" # for focal mean
  if (method == "max") {
    kernel[which(kernel != 0)] <- 1
  }
  # Apply focal function and write output
  copt <- c("COMPRESS=ZSTD", "ZLEVEL=1", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  opts <- list(gdal = copt, NAflag = nodata)
  terra::focal(dtm, kernel, method, filename = outpath, overwrite = TRUE, wopt = opts)
}

arguments <- commandArgs(trailingOnly = TRUE)
smooth(arguments[1], arguments[2])
