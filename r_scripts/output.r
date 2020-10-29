
########################################################################
#
# Towards a workflow for operational mapping of Aedes aegypti at urban 
# scale based on remote sensing
#
# Script written by: Pablo F. Cuervo and Veronica Andreo
#
########################################################################


#
# NOTES:
#
# This R script documents the steps followed to plot output distibution 
# maps for *Aedes aegypti* in urban areas.
#


#
# Set environment and load packages
#


# Set dependencies
library(checkpoint)
checkpoint("2020-07-01")

# Install required R packages (from CRAN)
install.packages("tmap", repos = "https://cran.rstudio.com")
install.packages("rgdal", repos = "https://cran.rstudio.com")
install.packages("raster", repos = "https://cran.rstudio.com")
install.packages("ymap", repos = "https://cran.rstudio.com")

# Load R packages
library(rgdal)
library(raster)
library(tmap)


#
# Set the working directory
#


# Set the working space
setwd("your_directory\\Project_Urban_Aedes")

# Load a polygon of the neighbourhoods in city of Cordoba
nhbrhd <- rgdal::readOGR("\\variables","barrios_cba")
nhbrhd <- sp::spTransform(nhbrhd, sp::CRS("+init=epsg:32720"))
nhbrhd$cba <- as.factor("cdad_cba")

# Load average and uncertainty raster from Maxent output (Example for December 2017)
dec17_mean <- raster::raster(list.files("\\Final_Models\\M_6_F_l_predictors_EC", 
                            pattern = glob2rx("*Cordoba_avg.asc"), 
                            full.names = TRUE))
raster::crs(dec17_mean) <- raster::crs(nhbrhd)

dec17_sd <- raster::raster(list.files("\\Final_Models\\M_6_F_l_predictors_EC", 
                           pattern = glob2rx("*Cordoba_stddev.asc"), 
                           full.names = TRUE))
raster::crs(dec17_sd) <- raster::crs(nhbrhd)


###
# Maps for Aedes aegypti 
###


#
# Average predicted probability and standard deviation
#


# average prediction
mean_plot <- 
    tmap::tm_shape(dec17_mean) + 
    tmap::tm_raster(palette = "YlOrRd", 
                    breaks = c(0,0.2,0.4,0.6,0.8,1)) +
    tmap::tm_shape(nhbrhd) + 
    tmap::tm_borders("gray20", lwd = .5) +
    tmap::tm_layout('December 2017 - Average Prediction', 
                    inner.margins=c(0,0,.1,0), 
                    title.size=.8, 
                    title.position = c('center', 'TOP')) + 
    tmap::tm_legend(show = FALSE)

# standard deviation
sd_plot <- 
    tmap::tm_shape(dec17_sd) + 
    tmap::tm_raster(palette = "YlGnBu", 
    breaks = c(0,0.05,0.1,0.15,0.2,0.3)) +
    tmap::tm_shape(nhbrhd) + 
    tmap::tm_borders("gray20", lwd = .5) +
    tmap::tm_layout('December 2017 - Standard Deviation', 
                    inner.margins=c(0,0,.1,0), 
                    title.size=.8, 
                    title.position = c('center', 'TOP')) +
    tmap::tm_legend(show = FALSE)

# map legends
mean_legend <- 
    tmap::tm_shape(dec17_mean) + 
    tmap::tm_raster(legend.is.portrait = FALSE, 
                    style="cont", palette = "YlOrRd", 
                    breaks = c(0,0.2,0.4,0.6,0.8,1)) +
    tmap::tm_shape(nhbrhd) + 
    tmap::tm_borders("gray20", lwd = .5) +
    tmap::tm_layout(legend.only = TRUE, 
                    legend.outside.position = c('center', 'TOP'), 
                    legend.title.color = 'white')
                    
sd_legend <- 
    tmap::tm_shape(dec17_sd) + 
    tmap::tm_raster(legend.is.portrait = FALSE, 
                    style="cont", palette = "YlGnBu", 
                    breaks = c(0,0.05,0.1,0.15,0.2,0.3)) +
    tmap::tm_shape(nhbrhd) + 
    tmap::tm_borders("gray20", lwd = .5) + 
    tmap::tm_layout(legend.only = TRUE, 
                    legend.outside.position = c('center', 'TOP'), 
                    legend.title.color = 'white')

# save plots
dev.new(height = 3.25, width = 5.5)
tmap::tmap_arrange(mean_plot, sd_plot, mean_legend, sd_legend, 
                   ncol=2, nrow=2, heights = c(0.75,0.25))


#
# Aggregated probabilities by neighborhood
#


# aggregate raster values in nhbrhd polygons
nhbrhd <- raster::extract(dec17_mean, nhbrhd, fun = mean, na.rm = TRUE)

# create figure
fig_neighb <- 
  tm_shape(nhbrhd) + 
  tm_fill(dec17_mean,
          palette = "plasma",
          n = 10, 
          colorNA = "white",
          legend.reverse = TRUE,
          title = "Probability") +
  tm_borders(lwd = 0.3) +
  tm_scale_bar(breaks = c(0, 3, 6), 
               text.size = 0.3,
               position = "right")

# save plot
tmap_save(fig_neighb, 
          filename= "fig_neighb.png", 
          width = 60,
          height = 45,
          units = "mm")


#
# Presence - Absence map
#


# threshold (obtained after the threshold dependent validation)
thrs <- 0.5

# apply threshold to average prediction
apply_thres <- function(x) { 
    ifelse(x <= thrs, 0,
    ifelse(x > thrs, 1, NA)) }

binary_map <- raster::calc(dec17_mean, fun=apply_thres)

# create figure
fig_binary_map <-  
    tmap::tm_shape(binary_map) + 
    tmap::tm_raster(style = "cat",
            palette = "-RdBu",
            n = 2,
            title = "Presence",
            labels = c("0", "1")) +
    tmap::tm_shape(nhbrhd) + 
    tmap::tm_borders("gray20", 
            lwd = .5) +
    tm_scale_bar(breaks = c(0, 3, 6), 
            text.size = 0.3,
            position = "right")

# save plot
tmap_save(fig_binary_map, 
          filename= "fig_binary_map.png", 
          width = 100,
          height = 80,
          units = "mm")


