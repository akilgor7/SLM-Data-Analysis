##
## Created 2020 07 30
## Last Edited 2022 03 17

library(mapview)
library(raster)
library(sp)
library(sdm)
library(usdm)
library(rgdal) # package will be phased out in 2023!
library(dismo)
library(maptools) # package will be phased out in 2023!
library(ggplot2)
library(bootStepAIC)
library(MASS)
library(GGally)

# dev.off()
installAll() #installs all associated with sdm package

setwd("~/Desktop/REPS:PuertoRico/LandslideAnalysis")

###### READING IN OCCURRENCE DATA #########
# pres = presence data
pres <- read.csv(file.path("2261-ii_point_data", "clip_1998_2006_int.csv"))
head(pres)

# add occurrence column
occurrence_pres <- c(rep(1, nrow(pres)))
pres <- cbind(occurrence_pres, pres)
# renaming occurrence column to have a standard for unified table with absence points
colnames(pres)[1] <- "occurrence"
head(pres)

## Add coordinate location & plotting points by location
coordinates(pres)=~lon+lat
plot(pres)
class(pres) #spatial points data frame

# create projection variable
CRS_points <- CRS("+proj=utm +zone=14 +datum=WGS84") # verify that this is correct projection
proj4string(pres) <- CRS_points
# proj4string(pres) <- projection(raster()) # try out from Marconi's lab

####### OBTAIN RASTER DATA #########
path <- system.file("crop_predictors", package="sdm")
lst <- list.files(path="crop_predictors", pattern = '.asc$', full.names=T)
lst

# matching extents
aspect <- raster(lst[1])
slope <- raster(lst[24])
elev <- raster(lst[22])
curv <- raster(lst[21])

bio1 <- raster(lst[2])

aspect <- setExtent(aspect, bio1, keepres=TRUE)
slope <- setExtent(slope, bio1, keepres=TRUE)
curv <- setExtent(curv, bio1, keepres=TRUE)
elev <- setExtent(elev, bio1, keepres=TRUE)

# combine raster layers
bio_lst <- lst[2:20]
bio_lst <- stack(bio_lst)

topo_lst <- c(aspect, slope, curv, elev)
  
pet <- raster(lst[23])

preds <- stack(bio_lst, pet, aspect, slope, curv, elev)

preds <- stack(preds, topo_lst) # this is repetative?
head(preds)
summary(preds)

###### RANDOM POINT SETUP ##############

# find boarders of 1991_2261-ii
max.lat <- max(pres$lat)
min.lat <- min(pres$lat)
max.lon <- max(pres$lon)
min.lon <- min(pres$lon)

#create extent
ext <- extent(min.lon, max.lon, min.lat, max.lat)

## create raster layer
mask_p <- pet

## check fit
plot(mask_p)
plot(ext, col='red', add=TRUE)

## create absence points
numPoints <- length(pres)
numPoints

set.seed(0)
backgr <- randomPoints(mask=mask_p, n=numPoints, p=pres, ext=ext, excludep=TRUE)
class(pres) #sp
class(backgr) #should be sp, same class is important

# add occurrence column
occurrence_abs <- c(rep(0, nrow(backgr)))
backgr <- cbind(occurrence_abs, backgr)
# renaming occurrence column
# standard column name with presence points
colnames(backgr)[1] <- "occurrence"
head(backgr)

# converting background points to data frame
backgr <- as.data.frame(backgr)
class(backgr)

# converting backgr to Spatial points data frame
backgr_coords <- backgr[,c(2,3)]
absence <- SpatialPointsDataFrame(coords = backgr_coords, data = backgr,
                proj4string = CRS_points)

# check absence points
plot(mask_p)
plot(ext, col='green', add=TRUE)
points(absence, cex=0.5, col='red')
points(pres, cex=0.5, col='blue')

##### COMBINE ABSENCE AND PRESENCE POINTS ######

# check projections
proj4string(pres) # getting a warning here
proj4string(absence) # same warning
# matching projections
pres <- spTransform(pres, proj4string(absence)) # same warning

# ERROR from lines 131/132: 
# Warning message:
# In proj4string(pres) :
#  CRS object has comment, which is lost in output; in tests, see
# https://cran.r-project.org/web/packages/sp/vignettes/CRS_warnings.html

# keeping occurrence column, removing everything else
pres_occurrence <- pres[,1]
abs_occurrence <- absence[,1]

# checking locations again
plot(pres_occurrence)
plot(abs_occurrence)

# merging presence and absence
sdm_occurrence <- rbind(pres_occurrence, abs_occurrence)
summary(sdm_occurrence)
head(sdm_occurrence)

# at this point check to validate:
#   - there should only be one column, occurrence
#   - mean value should be 0.5

### There are a few presence points with NA values for all raster info, what to do?
### should they be removed, and do they impact the sdm?

################################################################################
##### CHECKING FOR COLLINEARITY
################################################################################
v1 <- vifstep(preds)
v1

# keep what is not collinear
pcor <- exclude(preds, v1)
head(pcor)

### Why did we decide on v1 over v2? research this, make sure it is best method to deal 
### with collinearity paired with an sdm
################################################################################
##### MODEL FITTING 
################################################################################

# creating an sdm object
d <- sdmData(formula=occurrence~., train = sdm_occurrence, predictors= pcor)
class(d)

# using "d" to create a model with multiple methods
m1 <- sdm(occurrence~., data=d, methods=c("glm","rf","svm","brt","cart","mars"),test.percent=30, replication=c("sub"),n=3)
m1
# results from m1
# methods    :     AUC     |     COR     |     TSS     |     Deviance 
-------------------------------------------------------------------------
# glm        :     0.74    |     0.42    |     0.4     |     1.19     
# rf         :     0.92    |     0.75    |     0.72    |     0.71     
# svm        :     0.82    |     0.57    |     0.55    |     1.06     
# brt        :     0.82    |     0.54    |     0.55    |     1.14     
# cart       :     0.8     |     0.53    |     0.52    |     1.17     
# mars       :     0.79    |     0.5     |     0.49    |     1.19     

# view model in interface
gui(m1)

# Get model info
getModelInfo(m1)

# Get Roc curve
roc(m1)
roc(m1,smooth=T)

# Predictions
p1 <- predict(m1,pcor,mean=T) # creates a prediction for all 6 models
names(p1)
# PREDICTIONS FOR ALL MODEL METHODS USED: 
# [1] "sp_1.m_glm.re_subs"  "sp_1.m_rf.re_subs"   "sp_1.m_svm.re_subs"  "sp_1.m_brt.re_subs" 
# [5] "sp_1.m_cart.re_subs" "sp_1.m_mars.re_subs"
plot(p1$sp_1.m_glm.re_subs) # glm prediction, for example

################################################################################
##### Ensemble model
################################################################################

# Current distribution
# ensemble mixes all methods used to create m1

en=ensemble(m1,newdata=pcor,setting=list(method="weighted",stat="TSS",opt=2))
class(en)
plot(en,zlim=c(0,1))

# Ensemble for individual species or for specific models and models outputs
##Current
ev = getEvaluation(m1,stat=c("AUC","TSS","threshold"),opt=2)
mean(ev$threshold)

#### Download future scenario
biof=raster::getData("CMIP5",var="bio",res=0.5,lon=-66,lat=18,rcp=85,year=70,model="AC")
# add future data for Sierras de las Minas in place here
# under different scenarios?
names(preds)

names(biom)

names(biof)=names(bio)

date()
biof2=brick(biof$bio2_23,biof$bio3_23,biof$bio4_23,biof$bio8_23,biof$bio13_23,biof$bio14_23,biof$bio15_23,biof$bio18_23,biof$bio19_23)
date()
biof2

plot(biof)

plot(biof,xlim = c(-67.3,-65.5),ylim = c(17.8,18.6))
biof2
biom

###Predictions for future

p=predict(m1,biof,mean=T)


#### Ensemble

enF=ensemble(m1,newdata=biof2,setting=list(method="weighted",stat="TSS",opt=2))
plot(enF,xlim = c(-67.3,-65.5),ylim = c(17.8,18.6),col=tempcol(5),zlim=c(0,1))

### Comparing current and future
plot(stack(en,enf))






