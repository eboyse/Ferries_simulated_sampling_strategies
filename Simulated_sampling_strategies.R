#Code for 'Sampling from commercial vessel routes can capture marine biodiversity distributions effectively'

#download occurrence data from online repositories GBIF and OBIS 
library("rgbif")
library("robis")
library("dismo")
library("raster")
library("sp")
library("SSDM")

species.list <- read.csv("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/med.species.list.csv") #import list of species (latin names) of interest

species.list <- as.vector(species.list[5:14,1]) #create character vector of the latin species names

med.shape <- shapefile("C:/Users/betty/OneDrive/Documents/Leeds/Med Meta-analysis/med.shape/iho.shp") #shapefile to crop samples to area of interest

med.shape.EPSG3035 <- spTransform(med.shape, crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) #projection of interest

#function to download data from GBIF
download.from.gbif <- function(specieslist){
  for (i in 1:length(specieslist)){
    my.data <- occ_search(scientificName = as.character(specieslist[i]), return = "data", limit = 200000) ###download data, be careful of download limit
    my.data <- as.data.frame(my.data)
    my.data <- my.data[-c(which(is.na(my.data$decimalLatitude)),
                          which(is.na(my.data$decimalLongitude))), ] #removes data with NA for lon and lat
    my.data <- SpatialPointsDataFrame(data = my.data, coords = cbind(my.data$decimalLongitude, my.data$decimalLatitude),
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) #convert into spatial points object
    
    my.data <- spTransform(my.data, crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
    my.data <- crop(my.data, med.shape.EPSG3035)
    assign(paste0(sub(" ", ".", specieslist[i]), ".gbif"), my.data, env = globalenv())    
    write.csv(my.data, paste0("C:/Users/betty/OneDrive/Documents/Leeds/Med Meta-analysis/med.species.data.EPSG3035/", sub(" ", ".", specieslist[i]), ".gbif.csv")) #save .csv file with downloaded data 
    
    
    
  }
}


download.from.gbif(species.list)

#function to download data from OBIS
download.from.obis <- function(specieslist){
  for (i in 1:length(specieslist)){ 
    my.data <- occurrence(scientificname = as.character(specieslist[i]))
    my.data <- SpatialPointsDataFrame(data = my.data, coords = cbind(my.data$decimalLongitude, my.data$decimalLatitude),
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    my.data <- spTransform(my.data, crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
    my.data <- crop(my.data, med.shape.EPSG3035)
    assign(paste0(sub(" ", ".", specieslist[i]), ".obis"), my.data, env = globalenv())
    write.csv(my.data, paste0("C:/Users/betty/OneDrive/Documents/Leeds/Med Meta-analysis/med.species.data.EPSG3035/", sub(" ", ".", specieslist[i]), ".obis.csv"))
    
    
  }
}


download.from.obis(species.list)


#clean up data from GBIF and OBIS

#GBIF data 
setwd("C:/Users/betty/OneDrive/Documents/Leeds/Med Meta-analysis/med.species.data.EPSG3035")

species.list = list.files(pattern="*gbif.csv", full.names = T) #read in individual species .csv files and make one file for all GBIF entries
all.species <- read.csv(species.list[1])
all.species <- subset(all.species, select = c("species", "decimalLongitude", "decimalLatitude", "coords.x1", "coords.x2", "year", "month", "day", "datasetKey", "geodeticDatum", "basisOfRecord")) #keep data of interest
for (i in 2:length(species.list)) {
  work.species <- read.csv(species.list[i])
  work.species <- subset(work.species, select = c("species", "decimalLongitude", "decimalLatitude", "coords.x1", "coords.x2", "year", "month", "day", "datasetKey", "geodeticDatum", "basisOfRecord"))
  all.species <- rbind(all.species, work.species)
}

gbif.species <- subset(all.species, as.numeric(all.species$year) >= 2000) #remove records prior to 2000

#OBIS data - be aware that variable names are different between the datasets 
setwd("C:/Users/betty/OneDrive/Documents/Leeds/Med Meta-analysis/med.species.data.EPSG3035")

species.list = list.files(pattern="*obis.csv", full.names = T)
all.species <- read.csv(species.list[1])
all.species <- subset(all.species, select = c("scientificName", "decimalLongitude", "decimalLatitude", "coords.x1", "coords.x2", "date_year", "eventDate", "dataset_id", "geodeticDatum", "basisOfRecord"))
for (i in 2:length(species.list)) {
  work.species <- read.csv(species.list[i])
  work.species <- subset(work.species, select = c("scientificName", "decimalLongitude", "decimalLatitude", "coords.x1", "coords.x2", "date_year", "eventDate", "dataset_id", "geodeticDatum", "basisOfRecord"))
  
  all.species <- rbind(all.species, work.species)
}
obis.species <- subset(all.species, as.numeric(all.species$date_year) >= 2000)

#Accobams and Medlem datasets were manually downloaded 


###Stacked species distribution modelling - code for 'known' model 

###load environmental variables which were downloaded from Bio-Oracle or MarSpec 
bath <- raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/BO_bathymean_lonlat.tif")
tempmean <- raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/BO2_tempmean_ss_lonlat.tif")
chlomean <- raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/BO2_chlomean_ss_lonlat.tif")
temprange <- raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/BO2_temprange_ss_lonlat.tif")
bathslope <- raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/MS_biogeo06_bathy_slope_5m_lonlat.tif")
distshore <-raster("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Environmental_layers/env.var.med.WGS84/MS_biogeo05_dist_shore_5m_lonlat.tif")

#make raster stack - must all be in same resolution
env.pred <- stack(bath, tempmean, chlomean, temprange, bathslope, distshore)


##normalise env.var to between 0-1 
for (i in seq_len(length(env.pred@layers))){
  # For not categorical variable
  if (!env.pred[[i]]@data@isfactor) {
    env.pred[[i]] <- (env.pred[[i]] - env.pred[[i]]@data@min)/(env.pred[[i]]@data@max - env.pred[[i]]@data@min)
  }
}


####load occurrences
occurrences <- read.csv("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/med.species.data.EPSG6933.updated/rm3decplaces/ind_species_for_modelling/individual_species.csv") #code needs to include lon and lat coordinates

#spatially filter dataset using NND with distance 10 km which is often used for spatial filtering
NND.occurrences.list.10km <- thin(loc.data = occurrences, lat.col = "decimallatitude", long.col = "decimallongitude", 
                                  spec.col = "scientificnameaccepted", thin.par = 10, reps = 10000, 
                                  out.dir = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/spThin",
                                  locs.thinned.list.return = TRUE)

NND.occurrences <- as.data.frame(NND.occurrences.sp.10km)

#MAXENT modelling
#using default with NND occurrences 
MAX.default.NND <- modelling("MAXENT", Occurrences = NND.occurrences, Env = env.pred,
                             Xcol = "Longitude", Ycol = "Latitude",
                             cv = "holdout", cv.param = c(0.7, 3), bin.thres = "SES")


#Random forest modelling with 2 degree method
#using pseudo absences generated with 2 degree method (Barbet-Massin)
NND.circles <- circles(NND.occurrences.sp.10km, d = 222000)
NND.presence.area <- polygons(NND.circles) #2 degree area around presence points 

#randomly sample pseudo-absences outside of presence area
NND.pseudo.absence.area <- mask(env.pred[[1]], NND.presence.area, inverse = T) #mask present areas 

set.seed(1995)
NND.pseudo.absences <- randomPoints(NND.pseudo.absence.area, 20000) #randomly places points in areas not covered by presence points 
NND.pseudo.absences.sp <- SpatialPoints(coords = cbind(NND.pseudo.absences[,1], NND.pseudo.absences[,2]), crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) #make random points into spatial points object
NND.pseudo.absences.sp <- crop(NND.pseudo.absences.sp, med.shape)

NND.PA.sample.size <- length(NND.occurrences.sp.10km) #get number of presences 

NND.PA <- NND.pseudo.absences.sp[sample(1:length(NND.pseudo.absences.sp), NND.PA.sample.size),] #generate an equal number of pseudo absences as presences 

NND.pcol <- c(rep(1, nrow(NND.occurrences.sp.10km@coords)), rep(0, nrow(NND.PA@coords))) #make dataframe with column for presences (1) and absences (0) 

NND.occurrences_PA <- as.data.frame(rbind(NND.occurrences.sp.10km@coords, NND.PA@coords), header = T) #dataframe of presence and absence coordinates

NND.occurrences_PA$pcol <- NND.pcol #assign presence (1) and absence (0) column to the coordinates

#RF modelling  
RF.PA.NND <- modelling("RF", rf.args = list(ntree = 5000, nodesize = 5, do.trace = TRUE), Occurrences = NND.occurrences_PA, Env = env.pred,
                       Xcol = "coords.x1", Ycol = "coords.x2", Pcol = "pcol",
                       cv = "holdout", cv.param = c(0.7, 3), bin.thres = "SES")


#MARs modelling with 2 degree method
NND.MARS.circles <- circles(NND.occurrences.sp.10km, d = 222000)
NND.MARS.presence.area <- polygons(NND.MARS.circles)

NND.MARS.pseudo.absence.area <- mask(env.pred[[1]], NND.MARS.presence.area, inverse = T)

set.seed(1995)
NND.MARS.pseudo.absences <- randomPoints(NND.MARS.pseudo.absence.area, 20000) 
NND.MARS.pseudo.absences.sp <- SpatialPoints(coords = cbind(NND.MARS.pseudo.absences[,1], NND.MARS.pseudo.absences[,2]), crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
NND.MARS.pseudo.absences.sp <- crop(NND.MARS.pseudo.absences.sp, med.shape)

NND.MARS.PA.sample.size <- 1000 #1000 pseudo absences for all RF models 

NND.MARS.PA <- NND.MARS.pseudo.absences.sp[sample(1:length(NND.MARS.pseudo.absences.sp), NND.MARS.PA.sample.size),]

NND.MARS.pcol <- c(rep(1, nrow(NND.occurrences.sp.10km@coords)), rep(0, nrow(NND.MARS.PA@coords)))

NND.MARS.occurrences_PA <- as.data.frame(rbind(NND.occurrences.sp.10km@coords, NND.MARS.PA@coords), header = T)

NND.MARS.occurrences_PA$pcol <- NND.MARS.pcol

MARS.PA.NND <- modelling("MARS", Occurrences = NND.MARS.occurrences_PA, Env = env.pred,
                         Xcol = "coords.x1", Ycol = "coords.x2", Pcol = "pcol",
                         cv = "holdout", cv.param = c(0.7, 10), bin.thres = "SES")



###edited ensemble function from SSDM R package to calculate standard error for uncertainty instead of variance 
#############try and make new ensemble function but simply edit to SE
setGeneric("edited_ensemble", function(x, ..., name = NULL, ensemble.metric = c("AUC"),
                                       ensemble.thresh = c(0.75), weight = TRUE, thresh = 1001, uncertainty = TRUE, SDM.projections=FALSE, cores=0,
                                       verbose = TRUE, GUI = FALSE) {
  return(standardGeneric("edited_ensemble"))
})



###standard error function
standard_error <- function(x) sd(x) / sqrt(length(x))



setMethod("edited_ensemble", "Algorithm.SDM", function(x, ..., name = NULL, ensemble.metric = c("AUC"),
                                                       ensemble.thresh = c(0.75), weight = TRUE, thresh = 1001, uncertainty = TRUE, SDM.projections=FALSE, cores=0,
                                                       verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(name = name, ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, thresh = thresh, uncertainty = uncertainty, verbose = verbose,
             GUI = GUI)
  
  models <- list(x, ...)
  esdm <- Ensemble.SDM()
  sdms <- models
  
  # Algorithm ensemble model creation
  if (verbose) {
    cat("Creation of one ensemble niche model by algorithm...")
  }
  algo.ensemble <- list()
  while (length(models) > 0) {
    type.model <- list()
    type <- class(models[[1]])[[1]]
    rm <- {
    }
    for (i in seq_len(length(models))) {
      if (inherits(models[[i]], type)) {
        suppressWarnings({
          type.model[(length(type.model) + 1)] <- models[[i]]
        })
        rm <- c(rm, i)
      }
    }
    if (length(rm) > 0) {
      for (i in seq_len(length(rm))) {
        models[[rm[i]]] <- NULL
        rm <- rm - 1
      }
    }
    type.model["name"] <- name
    type.model[["ensemble.metric"]] <- ensemble.metric
    type.model[["ensemble.thresh"]] <- ensemble.thresh
    type.model["weight"] <- weight
    type.model["thresh"] <- thresh
    type.model["format"] <- TRUE
    type.model["verbose"] <- FALSE
    suppressWarnings({
      algo.ensemble[type] <- do.call(sum, type.model)
    })
  }
  if (verbose) {
    cat("   done. \n")
  }
  
  if (length(algo.ensemble) < 1) {
    if (verbose) {
      cat("No model were kept with this threshold, Null is returned. \n")
    }
    return(NULL)
  } else {
    # Extract indices of models that met the threshold
    ind <- 0
    selection.indices <- c()
    for (i in seq_len(length(sdms))) {
      if (all(sdms[[i]]@evaluation[,which(names(sdms[[i]]@evaluation) %in% ensemble.metric)] > ensemble.thresh)) {
        ind <- ind+1
        selection.indices[ind] <- i
      }
    }
    # Store individual models
    if(!SDM.projections){
      sdmsmin <- lapply(sdms[c(selection.indices)],function(x) {
        x@projection <- raster()
        x@binary <- raster()
        return(x)
      })
      esdm@sdms <- sdmsmin
    }
    else {
      esdm@sdms <- sdms[c(selection.indices)]
    }
    
    # Sum of algorithm ensemble
    if (verbose) {
      cat("Projection, evaluation and variable importance computing...")
    }
    algo.list <- list()
    for (i in seq_len(length(algo.ensemble))) {
      algo.list[[i]] <- algo.ensemble[[i]]
    }
    algo.list["name"] <- "sum"
    algo.list[["ensemble.metric"]] <- ensemble.metric
    algo.list[["ensemble.thresh"]] <- ensemble.thresh
    algo.list["weight"] <- weight
    algo.list["thresh"] <- thresh
    algo.list["format"] <- FALSE
    sum.algo.ensemble <- do.call(sum, algo.list)
    if (length(sum.algo.ensemble) < 1) {
      return(NULL)
    } else {
      
      # Name
      if (!is.null(name)) {
        name <- paste0(name, ".")
      } else {
        name <- "Specie."
      }
      esdm@name <- paste0(name, "Ensemble.SDM")
      
      # Projection
      esdm@projection <- sum.algo.ensemble@projection
      esdm@binary <- sum.algo.ensemble@binary
      if (verbose) {
        cat("   done \n")
      }
      
      # Data
      esdm@data <- algo.ensemble[[1]]@data
      if (length(algo.ensemble) > 1) {
        for (i in 2:length(algo.ensemble)) {
          esdm@data <- rbind(esdm@data, algo.ensemble[[i]]@data)
        }
      }
      
      # Evaluation
      if (verbose) {
        cat("Model evaluation...")
      }
      esdm@evaluation <- sum.algo.ensemble@evaluation
      if (verbose) {
        cat("   done \n")
      }
      
      # Axes evaluation
      if (verbose) {
        cat("Axes evaluation...")
      }
      esdm@variable.importance <- sum.algo.ensemble@variable.importance
      if (verbose) {
        cat("   done \n")
      }
      
      # Projections stack
      projections <- stack()
      for (i in seq_len(length(algo.ensemble))) {
        projections <- stack(projections, algo.ensemble[[i]]@projection)
        names(projections[[i]]) <- algo.ensemble[[i]]@name
      }
      
      # Algorithms Correlation
      if (!(uncertainty)) {
        if (verbose) {
          cat("Algorithm correlation computing is deactivated \n")
        }
      }
      if (uncertainty && length(projections@layers) > 1) {
        if (verbose) {
          cat("Algorithms correlation...")
        }
        esdm@algorithm.correlation <- as.data.frame(layerStats(projections,
                                                               "pearson", na.rm = TRUE)$`pearson correlation coefficient`)
        if (verbose) {
          cat("   done \n")
        }
      }
      
      # uncertainty map
      if (!(uncertainty)) {
        if (verbose) {
          cat("Uncertainty mapping is deactivated \n")
        }
      }
      if (uncertainty && length(projections@layers) > 1) {
        if (verbose) {
          cat("uncertainty mapping...")
        }
        if(cores>0){
          beginCluster(cores)
          esdm@uncertainty <- clusterR(projections, fun=function(x) calc(x, standard_error))
          endCluster()
        }
        else{
          esdm@uncertainty <- calc(projections, standard_error)
        }
        names(esdm@uncertainty) <- "uncertainty map"
        if (verbose) {
          cat("   done \n")
        }
      }
      
      # Algorithms Evaluation SDMs
      if (all(x@data$Presence %in% c(0, 1))) {
        if (verbose) {
          cat("Algorithms evaluation...")
        }
        esdm@algorithm.evaluation <- algo.ensemble[[1]]@evaluation
        row.names(esdm@algorithm.evaluation)[1] <- algo.ensemble[[1]]@name
        if (length(algo.ensemble) > 1) {
          for (i in 2:length(algo.ensemble)) {
            esdm@algorithm.evaluation <- rbind(esdm@algorithm.evaluation,
                                               algo.ensemble[[i]]@evaluation)
            row.names(esdm@algorithm.evaluation)[i] <- algo.ensemble[[i]]@name
          }
        }
        esdm@algorithm.evaluation$kept.model <- algo.ensemble[[1]]@parameters$kept.model
        if (length(algo.ensemble) > 1) {
          for (i in 2:length(algo.ensemble)) {
            esdm@algorithm.evaluation$kept.model[[i]] <- algo.ensemble[[i]]@parameters$kept.model
          }
        }
      } else {
        # MEMs
        warning("ALgorithms evaluation is not yet iplemented for MEMs !")
      }
      
      # Parameters
      esdm@parameters <- algo.ensemble[[1]]@parameters
      text.ensemble.metric <- character()
      text.ensemble.thresh <- character()
      for (i in seq_len(length(ensemble.metric))) {
        text.ensemble.metric <- paste0(text.ensemble.metric, ".",
                                       ensemble.metric[i])
        text.ensemble.thresh <- paste0(text.ensemble.thresh, "|",
                                       ensemble.thresh[i])
      }
      esdm@parameters$ensemble.metric <- text.ensemble.metric
      esdm@parameters$ensemble.thresh <- text.ensemble.thresh
      esdm@parameters$weight <- weight
    }
    
    if (verbose) {
      cat("   done \n")
    }
    
    return(esdm)
  }
})

#make ensemble model for each species 
Alopias_vulpinas <- edited_ensemble(MAX.default.NND, RF.PA.NND, MARS.PA.NND, SDM.projections = TRUE)
Alopias_vulpinas@name <- "Alopias_vulpinas" #name must be unique for each species for the stack function
save.esdm(Alopias_vulpinas, name = "Alopias_vulpinas", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/SSDMS_norm_env_var/Ensemble_sdms") #save to file 

#load individual species ensemble models
Alopias_vulpinas <- load_esdm(name = "Alopias_vulpinas", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/SSDMS_norm_env_var/Ensemble_sdms")

#make binary stacked species distribution model
bSSDM_stack <- stacking(Alopias_vulpinas, Balaenoptera_physalus, Carcharhinus_longimanus, Caretta_caretta, Cetorhinus_maximus, 
                        Conger_conger, Dasyatis_pastinaca, Delphinus_delphis, Dentex_dentex, Echelus_myrus, Echinorhinus_brucus,
                        Epinephelus_aeneus, Epinephelus_marginatus, Fistularia_commersonii, Globicephala_melas,
                        Grampus_griseus, Hexanchus_griseus, Isurus_oxyrinchus, Lophius_piscatorius, Merluccius_merluccius,
                        Mobula_mobular, Mola_mola, Molva_dypterygia, Muraena_helena, Myliobatis_aquila, Ophisurus_serpens,
                        Orcinus_orca, Physeter_macrocephalus, Pomatomus_saltatrix, Prionace_glauca, Raja_clavata,
                        Seriola_dumerili, Sphyraena_sphyraena, Sphyraena_viridensis, Squalus_acanthias, Stenella_coeruleoalba,
                        Thunnus_alalunga, Thunnus_thynnus, Torpedo_marmorata, Tursiops_truncatus, Xiphias_gladius,
                        Ziphius_cavirostris, Zu_cristatus, method = "bSSDM")

#save binary stacked species distribution model
save.stack(bSSDM_stack, name = "bSSDM_stack", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/SSDMS_norm_env_var")


###simulating different sampling strategies 

#load perfect knowledge SSDM (original SSDM made from online repository occurrence data using code above)
bSSDM <- load_stack(name = "bSSDM_stack", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/SSDMS_norm_env_var")

#get raster stack on binary species distribution rasters
pred <- stack(lapply(bSSDM@esdms, function(x){
  x@binary
}))

names(pred) <- unlist(lapply(strsplit(names(bSSDM@esdms), '.', fixed = TRUE), function(x){x[[1]]}))

pred <- projectRaster(pred, crs = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0", method = "ngb")


#load shapefile to be used as sampling frame 
all.ferries.sp <- shapefile("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Ferry_Routes_EPSG6933/all.ferries.shp")

#use sp function to similar different sampling strategies 
#need a spatial object as sampling frame, n = number of sampling points, type can be 'regular' or 'random' 
set.seed(3)
reg.samp.100 <- replicate(1000, spsample(all.ferries.sp, n=100, type='regular'), simplify = "matrix")

#select individual simulation, edit number between 1-40 to carry out models using different simulations 
regular.sample.sim1 <- reg.samp.100[[1]]

#get coordinates for sampling points 
reg.samp.sim1.coords <- as.data.frame(regular.sample.sim1@coords)

#extract occurrence data from original binary species distribution models for each sampling point
regular.sample.spec <- as.data.frame(raster::extract(pred, regular.sample.sim1))

species.names <- colnames(regular.sample.spec)

#loop to remake stacked species distribution models using the occurrence data from the sampling points
for (i in 1:length(regular.sample.spec)){
  species <- regular.sample.spec[i] #repeats for each species in stacked species distribution model
  species.coords <- cbind(species, reg.samp.sim1.coords)
  colnames(species.coords) <- cbind("species", "coords.x1", "coords.x2")
  species.pres <- species.coords[which(species.coords$species == "1"),]
  if (nrow(species.pres)>=20){ #minimum number of presences needed for each species for stacked species distribution model to be built
    species.sp <- SpatialPointsDataFrame(data = species.pres, coords = cbind(species.pres$coords.x1, species.pres$coords.x2),
                                         proj4string = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    species.WGS84 <- spTransform(species.sp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    species.final <- as.data.frame(species.WGS84@coords)
    species.final$species <- paste0(species.names[i])
    
    #Maxent 
    MAX.default <- modelling("MAXENT", Occurrences = species.final, Env = env.pred,
                             Xcol = "coords.x1", Ycol = "coords.x2",
                             cv = "holdout", cv.param = c(0.7, 3), bin.thresh = "SES")
    
    #Random forest
    circles <- circles(species.WGS84, d = 222000)
    presence.area <- polygons(circles)
    
    #randomly sample pseudo-absences outside of presence area
    pseudo.absence.area <- mask(env.pred[[1]], presence.area, inverse = T)
    
    set.seed(1995)
    pseudo.absences <- randomPoints(pseudo.absence.area, 20000) 
    pseudo.absences.sp <- SpatialPoints(coords = cbind(pseudo.absences[,1], pseudo.absences[,2]), crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    pseudo.absences.sp <- crop(pseudo.absences.sp, med.shape)
    
    PA.sample.size <- length(species.WGS84)
    
    PA <- pseudo.absences.sp[sample(1:length(pseudo.absences.sp), PA.sample.size),]
    
    pcol <- c(rep(1, nrow(species.WGS84)), rep(0, nrow(PA@coords)))
    
    occurrences_PA <- as.data.frame(rbind(species.WGS84@coords, PA@coords), header = T)
    
    occurrences_PA$pcol <- pcol
    
    RF.PA <- modelling("RF", rf.args = list(ntree = 5000, nodesize = 5, do.trace = TRUE), Occurrences = occurrences_PA, Env = env.pred,
                       Xcol = "coords.x1", Ycol = "coords.x2", Pcol = "pcol",
                       cv = "holdout", cv.param = c(0.7, 3), bin.thres = "SES")
    
    #MARS
    MARS.circles <- circles(species.WGS84, d = 222000)
    MARS.presence.area <- polygons(MARS.circles)
    
    MARS.pseudo.absence.area <- mask(env.pred[[1]], MARS.presence.area, inverse = T)
    
    set.seed(1995)
    MARS.pseudo.absences <- randomPoints(MARS.pseudo.absence.area, 20000) 
    MARS.pseudo.absences.sp <- SpatialPoints(coords = cbind(MARS.pseudo.absences[,1], MARS.pseudo.absences[,2]), crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    MARS.pseudo.absences.sp <- crop(MARS.pseudo.absences.sp, med.shape)
    
    MARS.PA.sample.size <- 1000
    
    MARS.PA <- MARS.pseudo.absences.sp[sample(1:length(MARS.pseudo.absences.sp), MARS.PA.sample.size),]
    
    MARS.pcol <- c(rep(1, nrow(species.WGS84@coords)), rep(0, nrow(MARS.PA@coords)))
    
    MARS.occurrences_PA <- as.data.frame(rbind(species.WGS84@coords, MARS.PA@coords), header = T)
    
    MARS.occurrences_PA$pcol <- MARS.pcol
    
    MARS.PA <- modelling("MARS", Occurrences = occurrences_PA, Env = env.pred,
                         Xcol = "coords.x1", Ycol = "coords.x2", Pcol = "pcol",
                         cv = "holdout", cv.param = c(0.7, 10), bin.thres = "SES")
    
    #make ensemble
    ensemble <- edited_ensemble(MAX.default, RF.PA, MARS.PA, SDM.projections = TRUE)
    ensemble@name <- paste0(species.names[i])
    save.esdm(ensemble, name = paste0(species.names[i]), path = paste0("C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Sample_SSDMs_NW_ferries/Reg_samp_100_occ20_1"))
    
  }
}

#open ESDMs into global environment 
file.names <- list.files(path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Sample_SSDMs_NW_ferries/Reg_samp_100_occ20_1")

for (i in 1:length(file.names)){
  esdm <- load_esdm(name = paste0(file.names[i]), path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Sample_SSDMs_NW_ferries/Reg_samp_100_occ20_1")
  assign(file.names[i], esdm)
}

#make binary stacked species distribution model
species.list <- lapply(file.names, as.name)

species.list <- paste0(lapply(file.names, as.name), collapse = ", ")

species.list #copy and paste output into stacking function 

reg.samp.100.SSDM <- stacking(Alopias_vulpinas, Balaenoptera_physalus, Carcharhinus_longimanus, Caretta_caretta, Cetorhinus_maximus, Conger_conger, Delphinus_delphis, Dentex_dentex, Globicephala_melas, Grampus_griseus, Hexanchus_griseus, Isurus_oxyrinchus, Lophius_piscatorius, Mobula_mobular, Mola_mola, Physeter_macrocephalus, Prionace_glauca, Stenella_coeruleoalba, Thunnus_thynnus, Torpedo_marmorata, Tursiops_truncatus, Xiphias_gladius, Ziphius_cavirostris, Zu_cristatus, 
                              method = "bSSDM")

save.stack(reg.samp.100.SSDM, name = "Reg_samp_100_occ20_1_SSDM", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Sample_SSDMs_NW_ferries")


###pearson correlation between original SSDM and remade SSDMs
original_bSSDM <- load_stack(name = "bSSDM_stack", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/SSDMS_norm_env_var") #perfect knowledge SSDM
remade_bSSDM <- load_stack(name = "Rand_samp_100_occ20_10_SSDM_random", path = "C:/Users/betty/OneDrive - University of Leeds/Leeds/Med Meta-analysis/Sample_SSDMS_15_routes") #load each remade SSDM individually

original_bSSDM_data <- getValues(original_bSSDM@diversity.map) #extract species richness from each raster cell
remade_bSSDM_data <- getValues(remade_bSSDM@diversity.map)
original_bSSDM_data <- na.omit(original_bSSDM_data)
remade_bSSDM_data <- na.omit(remade_bSSDM_data)
bSSDM_cor <- cor(original_bSSDM_data, remade_bSSDM_data) #pearson correlation coefficient 
bSSDM_cor

