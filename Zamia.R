# Owens Protocol for Workflows Project
# Zamia prasina

library(rangeModelMetadata)
library(occCite)
library(terra)
library(rangeBuilder)
library(geodata)
library(CoordinateCleaner)
library(voluModel)
library(ENMeval)
library(raster)
library(dplyr)
source("enm.maxnet.R") # From pending update to ENMeval, v2.0.5
#https://github.com/jamiemkass/ENMeval/blob/master/R/enm.maxnet.R

rmm <- rmmTemplate(family = c("base", "maxent"))
rmm$code$software$packages <- unique(rmmAutofillPackageCitation(rmm, (.packages()))$code$software$packages)
rmm$authorship$names <- 'Owens, Hannah L.'
rmm$authorship$contact <- 'hannah.owens@sund.ku.dk'
rmm$authorship$rmmName <- "Zamia prasina for Workflow Project"
rmm$authorship$license <- "CC"
rmm$authorship$contact <- "Hannah L. Owens, hannah.owens@sund.ku.dk"

rmm$studyObjective$purpose <- c("prediction","projection")
rmm$studyObjective$rangeType <- "potential"
rmm$studyObjective$invasion <- "native"
rmm$studyObjective$transfer <- "other"

# Get occurrences
gbifLoginHLO <- GBIFLoginManager()
occs <- occQuery("Zamia prasina", 
                 GBIFLogin = gbifLoginHLO, 
                 checkPreviousGBIFDownload = T, GBIFDownloadDirectory = "data/")
occRefs <- occCitation(occs)

rmm$data$occurrence$taxon <- occs@cleanedTaxonomy$`Best Match`
rmm$data$occurrence$dataType <- c("Presence only")
rmm$data$occurrence$sources <- sort(occRefs$occCitationResults$
                                      `Zamia prasina W.Bull`$Citation)

# Get environmental data ----
if(!file.exists("data/envtData/soil_world/phh2o_0-5cm_mean_30s.tif")){
  pH <- soil_world(var="phh2o", depth = 5, path = "data/envtData/")
} else pH <- rast("data/envtData/soil_world/phh2o_0-5cm_mean_30s.tif")

bioClim <- geodata::worldclim_global(var = "bio", res = 2.5, path = "data/envtData")
bioClim <- crop(bioClim,c(-100,-85,10,25))
bioClimVars <- c(1,6,12,14)
bioClim <- bioClim[[bioClimVars]]
if(!file.exists("data/envtData/soil_world/pHresampled.tif")){
  pHres <- resample(pH, bioClim, method = "bilinear", 
                    filename = "data/envtData/soil_world/pHresampled.tif")
} else {
  pHres <- rast("data/envtData/soil_world/pHresampled.tif")
}
preds <- c(bioClim,pHres)

rmm$data$environment$variableNames <- names(preds)
rmm$data$environment$sources <- c("@article{https://doi.org/10.1002/joc.5086,
  author = {Fick, Stephen E. and Hijmans, Robert J.},
  title = {WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas},
  journal = {International Journal of Climatology},
  volume = {37},
  number = {12},
  pages = {4302-4315},
  doi = {https://doi.org/10.1002/joc.5086},
  url = {https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/joc.5086},
  year = {2017}
}","@Article{soil-7-217-2021,
  AUTHOR = {Poggio, L. and de Sousa, L. M. and Batjes, N. H. and Heuvelink, G. B. M. and Kempen, B. and Ribeiro, E. and Rossiter, D.},
  TITLE = {SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty},
  JOURNAL = {SOIL},
  VOLUME = {7},
  YEAR = {2021},
  NUMBER = {1},
  PAGES = {217--240},
  URL = {https://soil.copernicus.org/articles/7/217/2021/},
  DOI = {10.5194/soil-7-217-2021}
}" 
)

# Future data
midCentury <- geodata::cmip6_world(model = "HadGEM3-GC31-LL", ssp = "585", time = "2041-2060", 
                                   var = "bioc", res = 2.5, path = "data/envtData/")[[bioClimVars]]

midCentury <- crop(midCentury, preds)
midCentury <- c(midCentury, preds$`phh2o_0-5cm`)
lateCentury <- geodata::cmip6_world(model = "HadGEM3-GC31-LL", ssp = "585", time = "2061-2080", 
                                    var = "bioc", res = 2.5, path = "data/envtData/")[[bioClimVars]]
lateCentury <- crop(lateCentury, preds)
lateCentury <- c(lateCentury, preds$`phh2o_0-5cm`)

# Enter metadata
rmm <- rmmAutofillEnvironment(rmm,preds,transfer=0) # for fitting environment
rmm <- rmmAutofillEnvironment(rmm,midCentury,transfer=1) # for transfer environment 1
rmm <- rmmAutofillEnvironment(rmm,lateCentury,transfer=2) # for transfer environment 2
rmm$data$environment$yearMin <- 1970
rmm$data$environment$yearMax <- 2000

# Clean occurrences and make training region
cleanOccs <- clean_coordinates(occs@occResults$`Zamia prasina W.Bull`$GBIF$OccurrenceTable,
                               tests = c("centroids", "equal", "gbif", 
                                         "institutions", "outliers", "zeros"),
                               lon = "longitude", lat = "latitude", species = "name", value = "clean", report = TRUE)

rmm$dataPrep$geographic$geographicalOutlierRemoval$rule <- "Remove if minimum distance to all other points is greater than 5 times interquantile range"
rmm$dataPrep$geographic$geographicalOutlierRemoval$notes <- "Default settings in coordinateCleaner::cleaned_coordinates()"
rmm$dataPrep$geographic$centroidRemoval$rule <- "Remove if <1000m from country centroid"
rmm$dataPrep$geographic$centroidRemoval$notes <- "Default settings in coordinateCleaner::cleaned_coordinates()"
rmm$dataPrep$biological$questionablePointRemoval$rule <- "latitude or longitude exactly 0; latitude or longitude an integer; latitude and longitude equal"
rmm$dataPrep$biological$questionablePointRemoval$notes <- "Default settings in coordinateCleaner::cleaned_coordinates()"

cleanOccs <- cleanOccs[any(as.numeric(cleanOccs$coordinateUncertaintyInMeters) > 5000,
                           is.na(cleanOccs$coordinateUncertaintyInMeters)),] # Remove high-uncertainty records
cleanOccs <- cleanOccs[cleanOccs$year > 1970,] # Remove particularly historical records
cleanOccs <- cleanOccs[cleanOccs$latitude > 15,] # Remove suspicious outlier after visual inspection
rmm$dataPrep$dataPrepNotes <- "Additional prep: removed occurrences with uncertainty greater than 5km (approximate resolution of environmental data; removed occurrences before 1970; removed occurrence at longitude -92.2732, latitude  14.9893 after manual inspection of remaining data"

rmm$data$occurrence$yearMin <- min(cleanOccs$year)
rmm$data$occurrence$yearMax <- max(cleanOccs$year)
rmm$data$occurrence$spatialAccuracy <- "Approximately 5000m uncertainty, maximum"

cleanOccs <- downsample(cleanOccs, preds[[1]]) # Thin to resolution of data
rmm$dataPrep$geographic$spatialThin$rule <- "Downsampled to resolution of training data (2.5 arcminutes)"
rmm$dataPrep$geographic$spatialThin$notes <- "Used voluModel::downsample()."

# Training region
trainingRegion <- vect(getDynamicAlphaHull(cleanOccs, coordHeaders = c("longitude", "latitude"), 
                                      fraction = 1, buff = 250000, verbose = T)[[1]])
rmm$data$environment$extentRule <- "Dynamic alpha hull fit around 100% of points (rangeBuilder::getDynamicAlphaHull()), buffered by 250km, unoccupied islands removed."
## Remove isolated spillover (e.g. nearby islands with no points)
occVect <- vect(cleanOccs, geom = c("longitude", "latitude"))
trainingRegion <- disagg(trainingRegion)
polysContainingPoints <- apply(relate(trainingRegion, occVect, "contains"),
                               MARGIN = 1, FUN = function(x) any(x))
trainingRegion <- aggregate(trainingRegion[polysContainingPoints])
trainPreds <- crop(preds, trainingRegion, mask = TRUE)
bg <- voluModel::mSampling2D(cleanOccs, trainPreds, mShp = trainingRegion)
set.seed(42)
bg <- bg[sample(1:nrow(bg), size = 10000, replace = F),]

# Make model ----
trainPreds <- stack(lapply(trainPreds, FUN = function(x) raster(x)))
model <- ENMevaluate(occs = cleanOccs, envs = trainPreds, bg = bg, 
                     tune.args = list(fc = c("L","LQ","LQP","Q", "QP","P"), rm = 1:3), 
                     partitions = "block", 
                     algorithm = "maxnet", doClamp = FALSE, rmm = rmm,
                     overlap = FALSE)
rmm$model$algorithms <- "maxnet"
rmm$model$speciesCount <- 1
rmm$model$algorithmCitation <- toBibtex(citation("maxnet"))

# Get fit statistics 
res <- eval.results(model)
opt.seq <- res %>% 
  filter(delta.AICc < 2) %>% 
  filter(auc.val.avg > (.9 * max(auc.val.avg))) %>% 
  filter(auc.diff.avg == min(auc.diff.avg))

# Select "best" model
modelOpt <- eval.models(model)[[opt.seq$tune.args]]
plot(modelOpt, type = "cloglog")
dev.off()
modProj <-  eval.predictions(model)[[opt.seq$tune.args]]
plot(modProj)
points(eval.bg(model), pch = 3, col = eval.bg.grp(model), cex = 0.5)
points(eval.occs(model), pch = 21, bg = eval.occs.grp(model))

# Project model ----
midCenturyTrain <- crop(midCentury, trainingRegion, mask = TRUE)
midCenturyTrain <- stack(lapply(midCenturyTrain, FUN = function(x) raster(x)))
names(midCenturyTrain) <- names(preds)
midCenturyProj <- maxnet.predict(m = modelOpt, envs = midCenturyTrain, 
                                 other.settings = model@other.settings)
midCenturySim <- similarity(ref = trainPreds, midCenturyTrain)$similarity_min
midCenturyProjNoMESS <- midCenturyProj * (midCenturySim >= 0)
writeRaster(rast(midCenturyProjNoMESS), "data/CycadMidCentury.tif")

lateCenturyTrain <- crop(lateCentury, trainingRegion, mask = TRUE)
lateCenturyTrain <- stack(lapply(lateCenturyTrain, FUN = function(x) raster(x)))
names(lateCenturyTrain) <- names(preds)
lateCenturyProj <- maxnet.predict(m = modelOpt, envs = lateCenturyTrain, 
                                  other.settings = model@other.settings)
lateCenturySim <- similarity(trainPreds, lateCenturyTrain)
lateCenturySim <- similarity(ref = trainPreds, lateCenturyTrain)$similarity_min
lateCenturyProjNoMESS <- lateCenturyProj * (lateCenturySim >= 0)
writeRaster(rast(midCenturyProjNoMESS), "data/CycadLateCentury.tif")