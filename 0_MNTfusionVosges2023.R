
################## MEAN MNT 
# Mean MNT on 25m, all files of 3 departements
############################################
#library(gdalUtils)
library(rgdal)
library(raster)
library(stringr)
library(terra)

datadir1 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D54\\D54\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D054_20211123"
MNT_files1 <- paste0(datadir1,"/",dir(datadir1, pattern = "\\.asc$"))

datadir2 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D67\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D067_20211123"
MNT_files2 <- paste0(datadir2,"/",dir(datadir2, pattern = "\\.asc$"))

datadir3 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D57\\D57\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D057_20211123"
MNT_files3 <- paste0(datadir3,"/",dir(datadir3, pattern = "\\.asc$"))

datadir4 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D68\\D68\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D068_20211123"
MNT_files4 <- paste0(datadir4,"/",dir(datadir4, pattern = "\\.asc$"))

datadir5 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D70\\D70\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D070_20210118"
MNT_files5 <- paste0(datadir5,"/",dir(datadir5, pattern = "\\.asc$"))

datadir6 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D88\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D088_20211123"
MNT_files6 <- paste0(datadir6,"/",dir(datadir6, pattern = "\\.asc$"))

datadir7 <- "Z:\\Data\\Vosges\\IGN\\RGEAlti\\D90\\D90\\RGEALTI\\1_DONNEES_LIVRAISON\\RGEALTI_MNT_1M_ASC_LAMB93_IGN69_D090_20210118"
MNT_files7 <- paste0(datadir7,"/",dir(datadir7, pattern = "\\.asc$"))



MNT_files <- c(MNT_files1,MNT_files2,MNT_files3,MNT_files4,MNT_files5,MNT_files6,MNT_files7)
#MNT_files <- MNT_files7
rm(MNT_files1,MNT_files2,MNT_files3,MNT_files4,MNT_files5,MNT_files6,MNT_files7)

#delete unuseful files not intersecting with studysite
library(sf)

# Load the study site polygon
study_site <- st_read('Z:\\Data\\Vosges\\D11limit35buffer.gpkg')  # Replace 'study_site.shp' with your polygon file
#plot(study_site$geom)

# List ASC files
asc_files <- MNT_files
rm(MNT_files)
# Loop through ASC files
study_extent <- extent(study_site)

for (asc_file in asc_files) {
  print(asc_file)
  # Load ASC file as an sf object (assuming it has spatial information)
  asc_data <- raster(asc_file)
  # plot(study_extent)
  # plot(asc_data,add=TRUE)
  crs(asc_data) <- "+init=epsg:2154"
  asc_extent <- extent(asc_data)

  # Check if it intersects with the study site polygon
  if (length(raster::intersect(asc_extent,study_site))<1) {
    # Delete the non-intersecting ASC file
    file.remove(asc_file)
    print("delete")
  }
}

asc_files <- paste0(datadir7,"/",dir(datadir7, pattern = "\\.asc$"))

for (asc_file in asc_files) {
  print(asc_file)
  # Load ASC file as an sf object (assuming it has spatial information)
  asc_data <- raster(asc_file)
  # plot(study_extent)
  # plot(asc_data,add=TRUE)
  crs(asc_data) <- "+init=epsg:2154"
  asc_extent <- ext(asc_data)
  asc_sf <- as.polygons(asc_extent)
  crs(asc_sf) <- "+init=epsg:2154"
  intersected <- st_intersection(st_as_sf(asc_sf),study_site)
  # Check if it intersects with the study site polygon
  if (nrow(intersected)==0) {
    # Delete the non-intersecting ASC file
    file.remove(asc_file)
    print("delete")
  }
}
# end intersection





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

ras_lst <- MNT_files
rm(MNT_files)


# # reduce number of significant numbers by rounding
# convertAllToRaster <- function(tif_file) {
#   r <- raster(tif_file)
#   ras2 <- round(r,2)
#   name <- str_sub(tif_file,-31,-5)
#   filename2 <- paste0("C:\\Users\\anouk.schleich\\Documents\\tempo\\",name,".tif")
#   # Give it a standardised name (r1,r2,r3, etc)
#   # Write the raster to file
#   raster::writeRaster(ras2, filename2,overwrite = T)
# }
# 
# 
# lapply(ras_lst, FUN = convertAllToRaster)
# 
# 
# convertAllToRaster <- function(tif_file) {
#   r <- raster(tif_file)
#   ras2 <- terra::aggregate(r, 2)
#   ras2 <- round(ras2,2)
#   name <- str_sub(tif_file,-31,-5)
#   filename2 <- paste0("C:\\Users\\anouk.schleich\\Documents\\tempo2\\",name,".tif")
#   # Give it a standardised name (r1,r2,r3, etc)
#   # Write the raster to file
#   raster::writeRaster(ras2, filename2,overwrite = T)
# }
# 
# lapply(ras_lst, FUN = convertAllToRaster)






# creation du virtual raster
t1 <-Sys.time()
terra::vrt(ras_lst, v, overwrite = T)
t2<-Sys.time()
ras <- rast(v)
t3<-Sys.time()

t3-t1



# preparation de la matrice pour la fonction focal

focal_weight_matrix  <- terra::focalMat(ras, 12.5, "circle")
#focal_weight_matrix[which(focal_weight_matrix!=0)] <-1    # si max

t4 <-Sys.time()

merged_raster_mean <- focal(ras, focal_weight_matrix, fun="sum")

#merged_raster_max <- focal(ras, focal_weight_matrix, fun="max") # si max
t5<-Sys.time()

t5-t4

merged_raster_mean <- raster(merged_raster_mean)

#merged_raster_max2 <- raster(merged_raster_max)

proj4string(merged_raster_mean) <- CRS("+init=epsg:2154")

#writeRaster(merged_raster_max, "X:\\Data\\Landes\\QuaspareDownloadIGN\\merged_raster_max_MNHC2022.tif" )

terra::writeRaster(merged_raster_mean, "Z:\\Data\\Vosges\\IGN\\RGEAlti\\RGEAlti_D11bufferMNTmean.tif")







