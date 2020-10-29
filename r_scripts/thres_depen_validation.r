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
# This R script documents the steps followed to perform threshold 
# dependent validation of distibution maps for *Aedes aegypti* in 
# urban areas.
#
# In this case we use positive and negative larvae records as 
# independent dataset to perform validation. 
#
# The script runs under R 3.6. Up to now, there's no update of
# SDMTools package to R 4.0.
#


#
# Set environment and load packages
#


# Set dependencies
library(checkpoint)
checkpoint("2020-07-01")

# Install required R packages (from CRAN)
install.packages("sp", repos = "https://cran.rstudio.com")
install.packages("rgdal", repos = "https://cran.rstudio.com")
install.packages("raster", repos = "https://cran.rstudio.com")
install.packages("SDMTools", repos = "https://cran.rstudio.com")

# Load R packages
library(sp)
library(rgdal)
library(raster)
library(SDMTools)


#
# Set the working directory and prepare data
#


# Set the working space
setwd("your_directory\\Project_Urban_Aedes")

# Load average raster from Maxent output (Example for December 2017)
dec17_mean <- raster::raster(list.files("\\Final_Models\\M_6_F_l_predictors_EC", 
                             pattern = glob2rx("*Cordoba_avg.asc"), 
                             full.names = TRUE))

# Read data on positive and negative records for larvae (subsequent month, i.e., Jan 2018)
neg_jan18 <- rgdal::readOGR("negatives_jan2018.shp")
pos_jan18 <- rgdal::readOGR("positives_jan2018.shp")

# create presence column with 0 and 1
neg_jan18$presence <- rep(0,length(neg_jan18$Name))
pos_jan18$presence <- rep(1, length(pos_jan18$Name))

# paste pos and neg together
larv_jan2018 <- rbind(pos_jan18,neg_jan18)

# reproject to utm20s (modify as needed according to the area)
larv_jan2018 <- sp::spTransform(larv_jan2018, sp::CRS("+init=epsg:32720"))

# extract predicted values for points of larval records
larv_jan2018 <- raster::extract(dec17_mean, larv_jan2018, 
                                method = 'bilinear', df = TRUE)


#
# Estimate threshold dependent measures
#


jan2018ot <- optim.thresh(larv_jan18$presence,larv_jan18$predicted)

accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$min.occurence.prediction)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$mean.occurence.prediction)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$`10.percent.omission`)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$`sensitivity=specificity`)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$`max.sensitivity+specificity`)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$max.prop.correct)
accuracy(larv_jan18$presence,larv_jan18$predicted,jan2018ot$min.ROC.plot.distance)

# estimate confusion matrix and other derived indices for each threshold
matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$min.occurence.prediction)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$mean.occurence.prediction)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$`10.percent.omission`)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$`sensitivity=specificity`)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$`max.sensitivity+specificity`)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$max.prop.correct)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))

matriz <- confusion.matrix(larv_jan18$presence,larv_jan18$predicted,
                           jan2018ot$min.ROC.plot.distance)
(FPR <- (matriz[2] / (matriz[1] + matriz[2]))) 
(FNR <- (matriz[3] / (matriz[3] + matriz[4])))
(Precision <- (matriz[4] / (matriz[2] + matriz[4])))
(Recall <- (matriz[4] / (matriz[3] + matriz[4])))
(NPV <- (matriz[1] / (matriz[1] + matriz[3])))
(OA <- ((matriz[1] + matriz[4]) / (matriz[1] + matriz[2] + matriz[3] + matriz[4])))
