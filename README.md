# GeoGEDI
Improving GEDI footprint geolocation using a Digital Elevation Model

Input needed : 
- GEDI L2A footprints : shotnumber, beam, delta_time, elev_lowestmode, lat_lowestmode, lon_lowestmode
- Digital Elevation Model of ground elevations, transformed to 25m mean raster

Script : 
- GeoGEDIalgorithmParallel.R : run GeoGEDI algorithm, parallelised by orbit.
- GeoGEDIalgorithm.R (old version : not parallelised)

Corresponding paper :
A. Schleich, S. Durrieu, M. Soma and C. Vega, "Improving GEDI Footprint Geolocation Using a High-Resolution Digital Elevation Model," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 16, pp. 7718-7732, 2023, doi: 10.1109/JSTARS.2023.3298991.

