# GeoGEDI
Improving GEDI footprint geolocation using a Digital Elevation Model

#### Research paper
A. Schleich, S. Durrieu, M. Soma and C. Vega, "Improving GEDI Footprint Geolocation Using a High-Resolution Digital Elevation Model," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 16, pp. 7718-7732, 2023, doi: 10.1109/JSTARS.2023.3298991.

## Input data
- Digital Elevation Model of ground elevations
- GEDI L2A footprints : shotnumber, beam, delta_time, elev_lowestmode, lat_lowestmode, lon_lowestmode

## Scripts
### DTMmean.R
Apply focal mean to input DTM (preferably VRT file).

#### Usage :
```bash
# Create a VRT with the right CRS and nodata (EPSG:5698 for France mainland, 5699 for Corsica)
gdalbuildvrt -a_srs EPSG:5698 -srcnodata -99999 D001.vrt RGEALTI_*_D001_*/**/*.asc
# How to loop over all departments
for directory in RGE* ; do
    dpt=$(echo $directory | cut -d_ -f6)
    gdalbuildvrt $dpt.vrt $directory/**/*.asc -srcnodata -99999 -a_srs EPSG:5698
done

### Run process
# Single file
./DTMmean.R D001.vrt
# Sequential
./DTMmean.R D001.vrt D002.vrt [...] D0X.vrt
# Async using GNU parallel
ls *.vrt | parallel -j 8 ./DTMmean.R {}
```

This will produce one file per VRT with name ${VRT basename}_smooth.tif, tiled and compressed with ZSTD level 1 + floating point predictor.

### GeoGEDIalgorithmParallel.R : run GeoGEDI algorithm, parallelised by orbit.
### GeoGEDIalgorithm.R (old version : not parallelised)
