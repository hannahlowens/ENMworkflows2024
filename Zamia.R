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
library(dplyr)

rmm <- rmmTemplate(family = c("base", "maxent"))
rmm$code$software$packages <- unique(rmmAutofillPackageCitation(rmm, (.packages()))$code$software$packages)
rmm$authorship$names <- 'Owens, Hannah L.'
rmm$authorship$contact <- 'hannah.owens@sund.ku.dk'
rmm$authorship$rmmName <- "Zamia prasina for Workflow Project"
rmm$authorship$license <- "CC"
rmm$authorship$contact <- "Hannah L. Owens, hannah.owens@sund.ku.dk"
rmm$authorship$authorNotes <- "None"
rmm$authorship$miscNotes <- "For Adam Smith workflows project"
rmm$authorship$doi <- "None"
rmm$authorship$relatedReferences <- "None"

rmm$studyObjective$purpose <- c("prediction","projection")
rmm$studyObjective$rangeType <- "potential"
rmm$studyObjective$invasion <- "native"
rmm$studyObjective$transfer <- "other;"

# Get occurrences
gbifLoginHLO <- GBIFLoginManager()
occs <- occQuery("Zamia prasina", 
                 GBIFLogin = gbifLoginHLO, 
                 checkPreviousGBIFDownload = T, GBIFDownloadDirectory = "data/")
occRefs <- occCitation(occs)

rmm$data$occurrence$taxon <- occs@cleanedTaxonomy$`Best Match`
rmm$data$occurrence$dataType <- c("presence only")
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
names(preds)[[5]] <- "surface_pH"
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
write.csv(cleanOccs, "data/CycadCleanOccs.csv", row.names = FALSE)

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
model <- ENMevaluate(occs = cleanOccs, envs = trainPreds, bg = bg, 
                     tune.args = list(fc = c("L","LQ","LQP","Q", "QP","P"), rm = 1:3), 
                     partitions = "block", 
                     algorithm = "maxnet", doClamp = FALSE, rmm = rmm,
                     overlap = FALSE)
rmm$model$algorithms <- "maxnet"
rmm$model$speciesCount <- 1
rmm$model$algorithmCitation <- toBibtex(citation("maxnet"))
rmm$model$partition$partitionRule <- "Spatial blocks defined by k means clustering, k = 4."
rmm$model$references <- c("@article{https://doi.org/10.1111/2041-210X.13628,
  author = {Kass, Jamie M. and Muscarella, Robert and Galante, Peter J. and Bohl, Corentin L. and Pinilla-Buitrago, Gonzalo E. and Boria, Robert A. and Soley-Guardia, Mariano and Anderson, Robert P.},
  title = {ENMeval 2.0: Redesigned for customizable and reproducible modeling of species’ niches and distributions},
  journal = {Methods in Ecology and Evolution},
  volume = {12},
  number = {9},
  pages = {1602-1608},
  doi = {https://doi.org/10.1111/2041-210X.13628},
  url = {https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13628},
  eprint = {https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210X.13628},
  year = {2021}
}")

# Get fit statistics and select "best" model
res <- eval.results(model)
opt.seq <- res %>% 
  filter(delta.AICc < 2) %>% 
  filter(auc.train > (.75 * max(auc.train))) %>% 
  filter(auc.diff.avg == min(auc.diff.avg))
write.csv(opt.seq, "data/LeopardCatFinalModelStats.csv", row.names = FALSE)
rmm$model$selectionRules <- "Successive filtering: delta AIC less then 2, models with AUC training scores at least 75% of maximum AUC training score, remaining model with the lowest AUCdiff"
rmm$assessment$trainingDataStats$AUC <- opt.seq$auc.train
rmm$assessment$testingDataStats$AUCDiff <- opt.seq$auc.diff.avg
rmm$assessment$trainingDataStats$AIC <- opt.seq$auc.diff.avg
rmm$model$finalModelSettings <- paste("Feature classes: ", opt.seq$fc[[1]], 
                                      "; Regularization multiplier: ", opt.seq$rm[[1]])

modelOpt <- eval.models(model)[[opt.seq$tune.args]]
plot(modelOpt, type = "cloglog")
rmm$assessment$notes <- "After examining response curves of preliminary plots, some environmental variables were removed due to lack of informativeness (i.e. they were nearly horizontal); statistics reported are for final model, AIC score reported is actually delta AICc."
dev.off()
modProj <-  eval.predictions(model)
modProj <- modProj[[names(modProj)== opt.seq$tune.args]]
plot(modProj)
points(eval.occs(model), pch = 21, bg = eval.occs.grp(model))
writeRaster(modProj, "data/CycadPresent.tif", overwrite = TRUE)

# Project model ----
rmm$prediction$extrapolation <- "extrapolate function"
rmm$prediction$transfer$notes <- "Climate model: HadGEM3-GC31-LL; mid-century: 2041-2060; late century: 2061-2080); Soil surface pH at future time points assumed to be the same as present."
midCenturyTrain <- crop(midCentury, trainingRegion, mask = TRUE)
midCenturyTrain <- c(midCenturyTrain, trainPreds$surface_pH)
names(midCenturyTrain) <- names(preds)
midCenturyProj <- maxnet.predictRaster(m = modelOpt, envs = midCenturyTrain, 
                                       other.settings = model@other.settings)
crs(midCenturyProj) <- crs(trainPreds)
midCenturySim <- similarity(ref = trainPreds, midCenturyTrain)$similarity_min
midCenturySim <- clamp(midCenturySim, lower=0, upper=Inf, values=FALSE)
midCenturyProjNoMESS <- mask(midCenturyProj, mask = midCenturySim)
writeRaster(midCenturyProjNoMESS, "data/CycadMidCentury.tif", overwrite = TRUE)

lateCenturyTrain <- crop(lateCentury, trainingRegion, mask = TRUE)
lateCenturyTrain <- c(lateCenturyTrain, trainPreds$surface_pH)
names(lateCenturyTrain) <- names(preds)
lateCenturyProj <- maxnet.predictRaster(m = modelOpt, envs = lateCenturyTrain, 
                                        other.settings = model@other.settings)
crs(lateCenturyProj) <- crs(trainPreds)
lateCenturySim <- similarity(ref = trainPreds, lateCenturyTrain)$similarity_min
lateCenturySim <- clamp(lateCenturySim, lower=0, upper=Inf, values=FALSE)
lateCenturyProjNoMESS <- mask(lateCenturyProj, mask = lateCenturySim)
writeRaster(lateCenturyProjNoMESS, "data/CycadLateCentury.tif", overwrite = TRUE)

# Tying a bow on the rmm
rmm$code$software$platform <- "R"
rmm$code$fullCodeLink <- "https://github.com/hannahlowens/ENMworkflows2024/blob/main/Zamia.R"

rmmCheckFinalize(rmm)
cleanRmm <- rmmCleanNULLs(rmm)

# Broken
rmmToCSV(cleanRmm, "CycadMetadata.csv")
