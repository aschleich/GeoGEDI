# GeoGEDI

Improving GEDI footprint geolocation using a Digital Elevation Model

## Research paper

A. Schleich, S. Durrieu, M. Soma and C. Vega, "Improving GEDI Footprint Geolocation Using a High-Resolution Digital Elevation Model," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 16, pp. 7718-7732, 2023, doi: 10.1109/JSTARS.2023.3298991.

## Requirements

You may use any OS with R 4.* and geospatial libraries installed (terra).  
A Linux OS is best recommended for the download and preprocessing steps.  
To do it using Docker, you can build an image as follow :  

```bash
docker build -t geogedi .
```

Default user in container will be geogedi (uid=1001). You can modify the "useradd" command in Dockerfile to match your user uid and gid.  

## Input data

- Digital Terrain Model
- GEDI L2A data: shotnumber, beam, delta_time, elev_lowestmode, lat_lowestmode, lon_lowestmode

## Scripts

### DTMmean.R

Apply focal mean to input DTM.

#### Input data

It is advised to use a VRT file instead of a huge GeoTIFF.  
Some preliminary steps may be required to obtain a single file input DTM :  

```bash
# Create a VRT with the right CRS and nodata (EPSG:5698 for France mainland, 5699 for Corsica)
gdalbuildvrt -a_srs EPSG:5698 -srcnodata -99999 D001.vrt RGEALTI_*_D001_*/**/*.asc
# How to loop over all departments
for directory in RGE* ; do
    dpt=$(echo $directory | cut -d_ -f6)
    gdalbuildvrt $dpt.vrt $directory/**/*.asc -srcnodata -99999 -a_srs EPSG:5698
done
```

#### Usage

The script takes 2 arguments :

- Input DTM (VRT or TIFF file)
- Output TIFF file path

```bash
# Should improve VRT read speed
export GDAL_DISABLE_READDIR_ON_OPEN=TRUE
# Simple
./DTMmean.R D0XX.vrt D0XX_smoothed.tif
# Using GNU parallel
bash -c "ls /data/*.vrt | parallel -j4 ./DTMmean.R {} outputs/{/.}.tif"
# Using a bash for loop in docker
docker run -e GDAL_DISABLE_READDIR_ON_OPEN=TRUE -v /data:/data geogedi bash -c "for f in /data/*.vrt ; do ./DTMmean.R $f /outputs/{f%%.vrt}.tif"
```

Output is a GeoTIFF, tiled and compressed with ZSTD (level 1 + floating point predictor).

### ExtractH5data.R

Preprocess a directory of HDF5 files and store usefull variables to a Rdata file.  

#### Input data

Input HDF5 files can be downloaded from EarthData website, using the GEDI02_A_V002 dataset located in EarthDataCloud.  
A bounding box ROI must be provided as your high res polygon would be transformed to a very simplifed shape.  
Spatial subset can be enabled in the processing settings, after products selection is done.  
Finer spatial subset will be performed using the R script.  

When the order is complete, you may download a list of URLs using the [EarthData dashboard](https://search.earthdata.nasa.gov/downloads).  
You'll need to setup a `.netrc` file following this [tutorial](https://harmony.earthdata.nasa.gov/docs#getting-started).  

The best way to download the products located in the Harmony endpoint is with curl. It is a fast server so you can do several downloads in parallel. Avoid using parallel while writing directly to a network storage, simply download to a fast storage and move it later.  

H5 files are uncompressed so it is advised to do it on-the-fly before writing to disk (~ 15x less disk space required).  

Use the following bash command to download and compress files, using curl and GNU parallel :  

```bash
# Install required CLI tools :
# sudo apt-get install curl parallel zip unzip
cat 3664487056-Harmony.txt | parallel -j 8 "curl -Lnbj --silent {} | zip -q > {/}.zip && printf \"@ -\n@={/}\n\" | zipnote -w {/}.zip"
```

The last part (printf and zipnote) is required in order to set the right filename within the zip file, since it is unknown ("-") because data is piped from stdin.  

#### Usage

You may then pre-process a group of files and output will be stored in a single Rdata file. H5 files are searched recursively in the input folder.  
The script takes 3 arguments :

- input directory of HDF5 files (unzipped)
- output Rdata file path
- (TODO) optional input polygon for spatial subset

```bash
cd /my/ouput/dir
./ExtractH5data.R <input H5 files directory> <output Rdata file> <shape-for-spatial-subset.gpkg>
```

### GeoGEDIalgorithmParallel.R

Run GeoGEDI algorithm, parallelized by orbit.  

### GeoGEDIalgorithm.R

Legacy script (to be removed).
