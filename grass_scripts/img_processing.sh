#!/usr/bin/env bash

########################################################################
#
# Towards a workflow for operational mapping of Aedes aegypti at urban 
# scale based on remote sensing
#
# Script written by: Veronica Andreo
#
########################################################################

#
# NOTES:
#
# This script assumes the user has GRASS GIS > 7.8 with all dependencies
# and extensions installed. See: https://grass.osgeo.org/download/
#
# Extensions needed for this script are:
#       i.wi
#       r.texture.tiled
#       r.diversity
#
# It also asumes that the grass database is set with locations and 
# mapsets created accordingly.
#
# The only step that requires intervention (so far) is the atmospheric 
# correction of SPOT bands. 
# See: https://grass.osgeo.org/grass78/manuals/i.atcorr.html
# for further details. Example of needed files are included in the repo.
#
# TODO:
# Implement this script with Sentinel 2 L2A data, to avoid the manual
# selection of dates and skip the atmospheric correction step.
#
# RUN this script as:
# `grass78 grassdata/locationname/mapsetname --exec sh img_processing.sh` 
#


#
# Set relevant variables 
#

export GRASS_OVERWRITE=1

DATA_FOLDER=/home/username/data
DATA_FILE=IMG_SPOT6_MS_201803151357271_ORT_C0000000033590_R1C1.TIF
RES=6

DATE=20180315

BLUE=SPOT_MS_${DATE}.B0
GREEN=SPOT_MS_${DATE}.B1
RED=SPOT_MS_${DATE}.B2
NIR=SPOT_MS_${DATE}.B3

# read from metadata of each scene
BLUE_GAIN=7.67 
GREEN_GAIN=9.25
RED_GAIN=10.34
NIR_GAIN=13.88
BIAS=0.0 # for all bands

TEXT_SIZE=33
PROCESSES=6

CLASS_NUM=15


#
# Import SPOT MS bands
#


# Band order in SPOT TIFF file:
# b2: red
# b1: green
# b0: blue
# b3: nir
# Values are ND (12 bits)
# GAIN and BIAS in metadata

r.import \
  resolution=value \
  resolution_value=${RES} \
  input=${DATA_FOLDER}/${DATA_FILE} \
  output=SPOT_MS_${DATE}


#
# Rename bands to avoid confusion (optional)
#


# band order in MS tif is: B2 (R), B1 (G), B0 (B), B3 (NIR)
g.rename raster=SPOT_MS_${DATE}.1,SPOT_MS_${DATE}.B2
g.rename raster=SPOT_MS_${DATE}.2,SPOT_MS_${DATE}.B1
g.rename raster=SPOT_MS_${DATE}.3,SPOT_MS_${DATE}.B0
g.rename raster=SPOT_MS_${DATE}.4,SPOT_MS_${DATE}.B3


#
# Set computational region to match imported data
#


g.region -p raster=${RED}


#
# Convert to TOA radiance: (DNb / GAINb) + BIASb
#


r.mapcalc expression="${BLUE}_rad = ${BLUE} / ${BLUE_GAIN}"
r.mapcalc expression="${GREEN}_rad = ${GREEN} / ${GREEN_GAIN}"
r.mapcalc expression="${RED}_rad = ${RED} / ${RED_GAIN}"
r.mapcalc expression="${NIR}_rad = ${NIR} / ${NIR_GAIN}"


#
# Atmospheric correction per band
#


for map in `g.list rast pat=*B?_rad` ; do 

  # get input range per band 
  eval `r.info -r $map` 
  
  # get band number
  band=`echo $map | cut -d. -f2 | cut -c 2`
  
  # atmospheric correction (parameter files must be ready)
  i.atcorr \
    input=$map \
    range=$min,$max \
    parameters=b${band}_spot_i.atcorr_param.txt \
    elevation=srtm_30m \
    output=$map_corr \
    rescale=0,1

done


#
# NDVI
#


i.vi \
  red=${RED}_rad_corr \
  nir=${NIR}_rad_corr \
  viname=ndvi \
  output=SPOT_MS_${DATE}_NDVI


#
# NDWI (mcfeeters)
#


i.wi \
  green=${GREEN}_rad_corr \
  nir=${NIR}_rad_corr \
  winame=ndwi_mf \
  output=SPOT_MS_${DATE}_NDWI


#
# Texture (over NIR band)
#


# add-on r.texture.tiled (parallelized version of r.texture)
g.extension extension=r.texture.tiled

# estimate texture measures
r.texture.tiled \
  input=${NIR}_rad_corr \
  size=${SIZE} \
  method=entr \
  tile_width=1000 \
  tile_height=1000 \
  processes=${PROCESSES} \
  output=${NIR}_entropy_${SIZE}

r.texture.tiled \
  input=${NIR}_rad_corr \
  size=${SIZE} \
  method=contrast \
  tile_width=1000 \
  tile_height=1000 \
  processes=${PROCESSES} \
  output=${NIR}_contrast_${SIZE}

r.texture.tiled \
  input=${NIR}_rad_corr \
  size=${SIZE} \
  method=corr \
  tile_width=1000 \
  tile_height=1000 \
  processes=${PROCESSES} \
  output=${NIR}_correlation_${SIZE}


#
# Create group of bands
#


i.group group=spot_bands \
  subgroup=all \
  input=`g.list type=raster pattern=SPOT* mapset=. sep=,`
  

#
# Signature files (stats for the classes)
#


i.cluster group=spot_bands \
  subgroup=all \
  signaturefile=signat_spot \
  classes=${CLASS_NUM} \
  reportfile=rep_spot.txt \
  separation=0.1 \
  iterations=40 \
  sample=15,15 --o


#
# Unsupervised classification
#


i.maxlik \
  group=spot_bands \
  subgroup=all \
  signaturefile=signat_spot \
  output=class_spot_${DATE}_${CLASS_NUM}c \
  reject=reject_class_spot_${DATE}_${CLASS_NUM}c


#
# Distance to different classes
#


for CLASS in `seq 1 15` ; do

  # create binary maps for each class
  r.mapcalc \
    expression="class_${CLASS} = if(class_spot_15c == ${CLASS}, class_spot_15c, null())"
  
  # distance
  r.grow.distance \
    input=class_${CLASS} \
    distance=distance_class_${CLASS}

done


#
# Class diversity: Simpson & Shannon indices
#


# add-on r.diversity
g.extension extension=r.diversity

r.diversity \
  input=class_spot_${DATE}_${CLASS_NUM}c \
  prefix=class_spot_${DATE}_${CLASS_NUM}c_diversity \
  size=${SIZE} \
  method=simpson

r.diversity \
  input=class_spot_${DATE}_${CLASS_NUM}c \
  prefix=class_spot_${DATE}_${CLASS_NUM}c_diversity \
  size=${SIZE} \
  method=shannon


#
# Context variables: richness, mode, interspersion
#


r.neighbors \
  input=class_spot_${DATE}_${CLASS_NUM}c \
  method=diversity \
  size=${SIZE} \
  output=class_spot_${DATE}_${CLASS_NUM}c_rich_${SIZE}

r.neighbors \
  input=class_spot_${DATE}_${CLASS_NUM}c \
  method=mode \
  size=${SIZE} \
  output=class_spot_${DATE}_${CLASS_NUM}c_mode_${SIZE}

r.neighbors \
  input=class_spot_${DATE}_${CLASS_NUM}c \
  method=interspersion \
  size=${SIZE} \
  output=class_spot_${DATE}_${CLASS_NUM}c_intersp_${SIZE}


#
# Mean and Standard deviation of NDVI and NDWI in neighborhoods
#


r.neighbors \
  input=SPOT_MS_${DATE}_NDVI \
  method=average \
  size=${SIZE} \
  output=SPOT_MS_${DATE}_NDVI_average_${SIZE}

r.neighbors \
  input=SPOT_MS_${DATE}_NDVI \
  method=stddev \
  size=${SIZE} \
  output=SPOT_MS_${DATE}_NDVI_sd_${SIZE}

r.neighbors \
  input=SPOT_MS_${DATE}_NDWI \
  method=average \
  size=${SIZE} \
  output=SPOT_MS_${DATE}_NDWI_average_${SIZE}

r.neighbors \
  input=SPOT_MS_${DATE}_NDWI \
  method=stddev \
  size=${SIZE} \
  output=SPOT_MS_${DATE}_NDWI_sd_${SIZE}


#
# Convert vector data to raster
#


v.to.rast -d \
  input=Canal \
  output=canal \
  use=val \
  value=1

v.to.rast -d \
  input=CursoAgua \
  output=curso_agua \
  use=val \
  value=1

v.to.rast -d \
  input=LineaFerrea
  output=linea_ferrea \
  use=val \
  value=1
  

#
# Estimate distance to water and railroads
#


r.grow.distance \
  input=canal \
  distance=distance_canal \
  metric=euclidean

r.grow.distance \
  input=curso_agua \
  distance=distance_curso_agua \
  metric=euclidean

r.grow.distance \
  input=linea_ferrea \
  distance=distance_linea_ferrea \
  metric=euclidean


#
# Export all maps outside GRASS
#


# distance to classes
for map in `g.list type=raster pattern=distance_class*` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# texture measures: Contr, Entr, Corr
for map in `g.list type=raster patttern="*_{Contr,Corr,Entr}"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# context info
for map in `g.list type=raster pattern="*_{intersp,mode,rich}*"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# diversity
for map in `g.list type=raster pattern="*{shannon,simpson}*"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# indices
for map in `g.list type=raster pattern="*ND[VW]I_"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# bands
for map in `g.list type=raster pattern="*MS*_corr"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# unsupervised classifications with 15 classes
for map in `g.list type=raster pattern="class_spot*_15c"` ; do
    r.out.gdal \
      input=${map} \
      output=${map}.tif \
      format=GTiff
done

# distance to water and railroads
for i in "distance_canal" "distance_curso_agua" "distance_linea_ferrea" ; do
    r.out.gdal \
      input=${i} \
      output=${i}.tif \
      format=GTiff
done


#
# Remove all maps from the mapset, except for DEM, distance to water and railroads
#

g.remove type=raster pattern=*${DATE}* -f
g.remove type=raster pattern=*class_* -f
g.remove type=group name=spot_bands
