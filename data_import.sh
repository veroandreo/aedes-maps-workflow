#!/usr/bin/env bash

########################################################################
#
# Towards a workflow for operational mapping of Aedes aegypti at urban 
# scale based on remote sensing
#
# Andreo et al.
#
########################################################################

#
# NOTES:
#
# This script assumes the user has GRASS GIS > 7.8 with all dependencies
# and extensions installed. See: https://grass.osgeo.org/download/
#
# Extensions needed for this script are:
#       r.in.srtm.region
#
# It also asumes that the grass database is set with locations and 
# mapsets created accordingly.
#
# This script is to be run only once since it imports ancillary data.
#
# User needs to be registered at https://urs.earthdata.nasa.gov/users/new
# In the user profile, two specific applications must be approved in 
# "My application" tab:
#
#    "LP DAAC Data Pool" application, and
#    "Earthdata Search" application.
#
# RUN this script as:
# `grass78 grassdata/locationname/mapsetname --exec sh data_import.sh` 
#


#
# Import SRTM DEM needed for atmospheric correction
#


# Install extension
g.extension r.in.srtm.region

# Proceed with SRTM download and import
r.in.srtm.region \
  output=srtm_30m \
  memory=300 \
  username=yourusername \
  password=yourpassword \
  method=bilinear \
  resolution=30


#
# Import relevant vector data for the city/area under study
#


# import neighborhoods
v.in.ogr \
  input=barrios.kml \
  output=neighborhoods \
  min_area=0.0001 \
  snap=-1

# import water canal
v.in.ogr \
  input=Canal.shp \
  output=Canal 

# import water courses  
v.in.ogr \
  input=Curso_Agua.shp \
  output=CursoAgua 

# import railroads
v.in.ogr \
  input=LineaFerrea.shp \
  output=LineaFerrea 


#
# Import urban area and estimate mean elevation (needed for atmospheric correction files)
#


v.in.ogr \
  input=planta_urbana.shp \
  output="planta_urbana"
  
v.rast.stats \
  map=planta_urbana \
  raster=srtm_30m \
  column_prefix=dem \
  method=average

v.db.select \ 
  map=planta_urbana \
  where="FNA == 'Gran Cordoba'"



























