
################## MEAN MNT 
# Mean MNT on 25m, all files of 3 departements
############################################
#library(gdalUtils)
library(rgdal)
library(raster)
library(stringr)
library(terra)
library(sf)


write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
Sys.which("make")
install.packages("sp", type = "source")
library(jsonlite)
library(raster)
library(sp)
paste(RTOOLS40_HOME)

install.packages("sp",INSTALL_opts="--no-multiarch")

################################### version VRT ##########################################


library(terra)
library(Rcpp)

#write("TMP = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
rasterOptions(tmpdir='Z:/Rtemp_Anouk')
rasterOptions()

terraOptions(tempdir='Z:/Rtemp_Anouk')
terraOptions()


v <-   "Z:\\Data\\Vosges\\IGN\\RGEAlti\\vrt_RGEAlti_D11buffer.vrt" #"path/new.vrt" # name of virtual raster

datadir1 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D54\\D54\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D054_20211123"
ras_lst <- paste0(datadir1,"/",dir(datadir1, pattern = "\\.asc$"))



# creation du virtual raster
terra::vrt(ras_lst, v, overwrite = T)
ras <- rast(v)



# preparation de la matrice pour la fonction focal

focal_weight_matrix  <- terra::focalMat(ras, 12.5, "circle")
#focal_weight_matrix[which(focal_weight_matrix!=0)] <-1    # if max


merged_raster_mean <- focal(ras, focal_weight_matrix, fun="sum")

#merged_raster_max <- focal(ras, focal_weight_matrix, fun="max") # if max


merged_raster_mean <- raster(merged_raster_mean)

proj4string(merged_raster_mean) <- CRS("+init=epsg:2154")

#writeRaster(merged_raster_max, "X:\\Data\\Landes\\QuaspareDownloadIGN\\merged_raster_max_MNHC2022.tif" ) #if max

terra::writeRaster(merged_raster_mean, "Z:\\Data\\Vosges\\IGN\\RGEAlti\\RGEAlti_D11bufferMNTmean.tif")







