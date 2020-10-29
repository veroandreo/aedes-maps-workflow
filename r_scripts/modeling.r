
########################################################################
#
# Towards a workflow for operational mapping of Aedes aegypti at urban 
# scale based on remote sensing
#
# Script written by: Pablo F. Cuervo
#
########################################################################


#
# NOTES:
#
# This R script documents the steps followed to develop distibution 
# maps for *Aedes aegypti* in urban areas.
#
# See https://github.com/marlonecobos/kuenm for details about 
# development of ecological niche models with Maxent using the "kuenm" 
# package.
#


#
# Set environment and load packages
#


# Set dependencies
library(checkpoint)
checkpoint("2020-07-01")

# Install required R packages (from CRAN)
install.packages("rgdal", repos = "https://cran.rstudio.com")
install.packages("sp", repos = "https://cran.rstudio.com")
install.packages("maptools", repos = "https://cran.rstudio.com")
install.packages("raster", repos = "https://cran.rstudio.com")
install.packages("snowfall", repos = "https://cran.rstudio.com")
install.packages("SDMtune", repos = "https://cran.rstudio.com")
install.packages("devtools", repos = "https://cran.rstudio.com")
install.packages("Epi", repos = "https://cran.rstudio.com")

# Install required R packages (from GitHub)
devtools::install_github("marlonecobos/kuenm")

# Load R packages
library(rgdal)
library(sp)
library(maptools)
library(raster)
library(snowfall)
library(SDMtune)
library(kuenm)
library(Epi)


#
# Set the working directory
#


# Set the working space
setwd("your_directory\\Project_Urban_Aedes")


#
# Define the working space
#


# In this case, the working space is limited to the city of Cordoba 
# (Argentina). You should define your own geographic space.

# Load a polygon of the neighbourhoods in city of Cordoba
nhbrhd <- rgdal::readOGR("\\variables","barrios_cba")
nhbrhd <- sp::spTransform(nhbrhd, sp::CRS("+init=epsg:32720"))
nhbrhd$cba <- as.factor("cdad_cba")

# Polygon of neighbourhoods
cba_city <- maptools::unionSpatialPolygons(nhbrhd, nhbrhd$cba)


#
# Prepare the data required
#


# Load oviposition data from local sanitary authorities
# In this case, entomological samples were collected and processed 
# by the Zoonosis Division of the Health Ministry of Córdoba province
# (Argentina).
Data <- read.table("Ovip_2017_2018.txt", sep = "\t", header = TRUE)

# Limit data to a week prior and a week after the Spot scene date
# (in this case, week 50 = Dec 2017)
Data <- Data[,c(1:4,14:16)] # Selection of weeks 49 t0 51

# Cumulative oviposition in the selected period (weeks 49-51)
Data$cum_ovip <- rowSums(Data[,c(5:7)], na.rm = TRUE)

# Define presence/absence of oviposition in each ovitrap
Data$pres_abs <- Data$cum_ovip
Data$pres_abs[Data$cum_ovip > 0] <- 1     # Oviposition greater than 0 equals "1" (presence)
Data$pres_abs[Data$cum_ovip == 0] <- 0    # No oviposition equals "0" (absence)
  
# Transform decimal coordinates to UTM to fit SPOT imagery
dec.coords <- sp::SpatialPoints(cbind(Data$Long, Data$Lat), proj4string = sp::CRS("+proj=longlat"))
UTM.coords <- sp::spTransform(dec.coords, sp::CRS("+init=epsg:32720"))
Data$Long.UTM <- (UTM.coords$coords.x1)
Data$Lat.UTM <- (UTM.coords$coords.x2)
Data <- Data[,c(1,2,4,3,10,11,5:9)]


#
# Define accessible area M considering the presence of eggs in ovitraps
#


# Define a buffer of 800 m radio (maximum dispersal distance) around 
# each presence record
pres <- Data[,c(5:6,11)]
M_buffer <- raster::buffer(sp::SpatialPoints(pres[pres$pres_abs==1, 1:2]), width = 800)


#
# Map of ocurrence localities and accessible area M (December 2017)
#


raster::plot(cba_city, axes = TRUE, col = "light yellow", xlab = "Longitude (m E)", ylab= "Latitude (m S)")
raster::plot(M_buffer, col = "lightblue", add = TRUE)
points(pres$Long.UTM,pres$Lat.UTM,col="orange",pch=20,cex=0.75)
points(pres$Long.UTM,pres$Lat.UTM,col="red",cex=0.75)


#
# Create datasets for model calibration
#


# Set presence data for calibration
aedes_ovip <- pres[pres$pres_abs == 1,]
aedes_ovip$species <- "Aedes aegypti"
aedes_ovip <- aedes_ovip[,c(4,1,2)]

names(aedes_ovip)[names(aedes_ovip)=="Long.UTM"] <- "longitude"
names(aedes_ovip)[names(aedes_ovip)=="Lat.UTM"] <- "latitude"

# Split presence data for training and testing 
# (with 25% of occurrences) in calibration process
set.seed(1)
fold.pres		<- kfold(aedes_ovip, k = 4)
traindata		<- aedes_ovip[fold.pres != 1, ]
testdata		<- aedes_ovip[fold.pres == 1, ]

# Save datasets for calibration
write.csv(ovip_aedes, "aedes_joint.csv", row.names = FALSE)
write.csv(traindata, "aedes_train.csv", row.names = FALSE)
write.csv(testdata, "aedes_test.csv", row.names = FALSE)


#
## Load environmental variables derived from satellite imagery (i.e. SPOT, Sentinel)
#


# List of environmental variables
files <- list.files("\\variables", pattern = "tif", full.names = TRUE)

# Stack of environmental variables as predictors
predictors <- stack(files)


#
# Save predictors to G and M
#

# Parallel processing
require(snow)
snowfall::sfInit(parallel = TRUE, cpus=parallel:::detectCores()-2) 
snowfall::sfLibrary(raster)

# Save environmental rasters to G as .asc (as required by the "kuenm" package)
dir.create('G_space')
dir.create('G_space/predictors')
nums <- names(predictors)
writeRaster(predictors, filename = "G_space/predictors/.asc", 
            format = "ascii", bylayer = TRUE,
            suffix = nums, prj = TRUE, overwrite = TRUE)

# Crop and mask layers to M
M_variables <- raster::mask(raster::crop(predictors, M_buffer), M_buffer)

# Save environmental rasters to M as .asc (as required by the "kuenm" package)
dir.create('M_space')
dir.create('M_space/predictors')
nums <- names(M_variables)
writeRaster(M_variables, filename = "M_space/predictors/.asc", 
            format = "ascii", bylayer = TRUE,
            suffix = nums, prj = TRUE, overwrite = TRUE)



###
# Workflow for the ecological niche modelling
####


#
# Calibration of preliminary models for variable selection
#


# Preparing variables to be used in arguments
kuenm::kuenm_start(file.name = "Preliminary_model")

occ_joint <- "aedes_joint.csv"
occ_tra <- "aedes_train.csv"
M_var_dir <- "M_space"
batch_cal <- "Preliminary_model"
out_dir <- "Preliminary_Model"
reg_mult <- c(seq(0.1,1,0.1), seq(2, 6, 1), 8, 10)	# regularization multiplier(s) to be evaluated
f_clas <- c('l','lq','lp','lqp','h','lh','lqh','lph','lqph') # feature class(es) to be evaluated
args <- NULL
maxent_path <- "your_path_to_the_maxent.jar_file"
wait <- FALSE
run <- TRUE

kuenm::kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, 
                 M.var.dir = M_var_dir, batch = batch_cal, 
                 out.dir = out_dir, reg.mult = reg_mult, 
                 f.clas = f_clas, args = args, 
                 maxent.path = maxent_path, wait = wait, run = run)


#
# Evaluation and selection of best preliminary models
#


# Parallel processing
require(snow)
snowfall::sfInit(parallel=TRUE, cpus=parallel:::detectCores()-3)

# Load the required packages inside the cluster
snowfall::sfLibrary(kuenm)

occ_test <- "aedes_test.csv"
out_eval <- "Calibration_results_preliminary"
out_dir <- "Preliminary_Model"
threshold <- 5
rand_percent <- 25
iterations <- 5000
kept <- FALSE
selection <- "OR_AICc"
paral_proc <- FALSE

kuenm::kuenm_ceval(path = out_dir, occ.joint = occ_joint, 
                   occ.tra = occ_tra, occ.test = occ_test, 
                   batch = batch_cal, out.eval = out_eval, 
                   threshold = threshold, rand.percent = rand_percent, 
                   iterations = iterations, kept = kept, 
                   selection = selection, parallel.proc = paral_proc)
snowfall::sfStop()


#
# Identification of the best preliminary model
#


# Identify selected Sets
selected_models <- read.csv("your_directory\\Project_Urban_Aedes\\Calibration_results\\best_candidate_models_OR_AICc.csv", sep = ",", header = TRUE)

Sets <- c()
for (i in 1:nrow(selected_models)) {
  Set <- selected_models[i, 1]
  Sets[i] <- gsub(".*(Set_)", "\\1", Set)
}
Sets


#
# Variable selection with the "SDMtune" package
#


# Select random points as background information
bg <- dismo::randomPoints(M_variables, 10000)

# Create SWD data
swd_data <- SDMtune::prepareSWD(species = "Aedes aegypti", p = aedes_ovip[,2:3], a = bg, env = M_variables)

# Train model with parameters from best preliminary model
Folds <- SDMtune::randomFolds(data = swd_data, k = 30, only_presence = TRUE, seed = 1)
model <- SDMtune::train(method = "Maxent", data = swd_data, folds = Folds, fc = "l", reg = 6, iter = 500)

# Remove highly correlated variables
Background.data <- SDMtune::prepareSWD(species = "Aedes aegypti", 
                                       a = bg, 
                                       env = M_variables)

selected_variables_model <- SDMtune::varSel(model, metric = "tss", 
                                            bg4cor = Background.data, 
                                            method = "spearman", 
                                            cor_th = 0.7, permut = 1)

# Remove variables with low percent contribution
reduced_variables_model <- SDMtune::reduceVar(selected_variables_model, 
                                              metric = "tss", th = 5, 
                                              permut = 1, use_pc = TRUE)

#
# The selected, uncorrelated and most relevant environmental variables 
# are manually copied to a new folder in the working directory named 
# "M_variables_selected"
#


#
# Calibrate final models with uncorrelated and most relevant environmental variables
#


# Prepare variables to be used in arguments
kuenm::kuenm_start(file.name = "Calibration_final_model")

occ_joint <- "aedes_joint.csv"
occ_tra <- "aedes_train.csv"
M_var_dir <- "M_variables_selected"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1,1,0.1), seq(2, 6, 1), 8, 10) # regularization multiplier(s) to be evaluated
f_clas <- c('l','lq','lp','lqp','h','lh','lqh','lph','lqph') # feature class(es) to be evaluated
args <- NULL
maxent_path <- "your_path_to_the_maxent.jar_file"
wait <- FALSE
run <- TRUE

kuenm::kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, 
                 M.var.dir = M_var_dir, batch = batch_cal, 
                 out.dir = out_dir, reg.mult = reg_mult, 
                 f.clas = f_clas, args = args, 
                 maxent.path = maxent_path, wait = wait, 
                 run = run)


#
# Evaluation and selection of best models 
#

# Parallel processing
require(snow)
snowfall::sfInit(parallel=TRUE, cpus=parallel:::detectCores()-3)

# Load the required packages inside the cluster
snowfall::sfLibrary(kuenm)

occ_test <- "aedes_test.csv"
out_eval <- "Calibration_results_final"
out_dir <- "Candidate_Models"
threshold <- 5
rand_percent <- 25
iterations <- 5000
kept <- FALSE
selection <- "OR_AICc"
paral_proc <- FALSE

kuenm::kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra,
                   occ.test = occ_test, batch = batch_cal, 
                   out.eval = out_eval, threshold = threshold, 
                   rand.percent = rand_percent, iterations = iterations, 
                   kept = kept, selection = selection, 
                   parallel.proc = paral_proc)

snowfall::sfStop()


#
# Create final models
#


batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 30
rep_type <- "Bootstrap"
jackknife <- TRUE
out_format <- "cloglog"
M_var_dir <- "M_variables_selected"
project <- TRUE
G_var_dir <- "G_variables"
ext_type <- "ext_clam"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- c("writebackgroundpredictions = true", "writeplotdata = true") 

kuenm::kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, 
                 out.eval = out_eval, batch = batch_fin, rep.n = rep_n, 
                 rep.type = rep_type, jackknife = jackknife, 
                 out.dir = mod_dir, out.format = out_format, 
                 project = project, G.var.dir = G_var_dir, 
                 ext.type = ext_type, write.mess = write_mess, 
                 write.clamp = write_clamp, maxent.path = maxent_path, 
                 args = args, wait = wait1, run = run1)


#
# Calculating AUC, Specificity and Sensitivity from Maxent outputs
#


outdir <- '\\Final_Models\\M_6_F_l_predictors_EC'

P <- list.files(outdir, 'samplePredictions\\.csv$', full.names=TRUE)
Presence <- list()
for(x in 1:length(P)) {
  d1 <- subset(read.csv(P[x]), Test.or.train == 'train')[,6]
  Presence[[x]] <- d1
  }
Presence <- unlist(Presence, use.names=F)

Bg <- list.files(outdir, 'backgroundPredictions\\.csv$', full.names=TRUE)
Background <- list()
for(x in 1:length(Bg)) {
  d0 <- subset(read.csv(Bg[x]))[,5]
  Background[[x]] <- d0
  }
Background <- unlist(Background, use.names=F)

# ROC
Epi::ROC(c(Presence, Background), 
         rep(1:0, c(length(Presence), length(Background))), 
         PS=T, 
         plot="ROC",
         MX=F, 
         MI=F)
