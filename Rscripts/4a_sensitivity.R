################################################################################################################################
#   Performs sensitivity analysis of time to hidraulic failure and time to stomatal closure in ", model, "-medfate
################################################################################################################################
library(medfate)
library(sensitivity)
library(ggplot2)
library(doParallel)
library(boot)

# remotes::install_github("emf-creaf/medfate", ref = "devel")

# Sensitivity parameters --------------------------------------------------
n <- 10 #1000 # Number of rows (combinations) in the parameter matrices
nboot <- 0 #10 # Number of bootstrap samples
ncores <- 6 #20 # Number of cores
model <- "Sperry"

# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270
pue_slope <- 0
pue_aspect <- 0


# Common soil -------------------------------------------------------------
common_soil <- function(){
  fb_soil_df <- defaultSoilParams(n = 3)
  fb_soil_df$widths <- c(200, 800, 3000)
  fb_soil_df$rfc <- c(30, 50, 94.5)
  fb_soil_df$clay <- rep(39, 3)
  fb_soil_df$sand <- rep(26, 3)
  fb_soil_df$om <- c(6,3,1)
  fb_soil_df$bd <-rep(1.45, 3)
  
  fb_soil <- soil(fb_soil_df)
  fb_soil$VG_theta_sat <- rep(0.5, 3)
  fb_soil$VG_theta_res <- rep(0.098, 3)
  fb_soil$VG_n <- rep(1.55, 3)
  # 1 cm = 0.00009804139432 MPa
  fb_soil$VG_alpha <- rep(0.0005/0.00009804139432, 3)
  fb_soil$Ksat <- rep(655.2934*1.69, 3)
  
  return(fb_soil)
}
soil <- common_soil()


# Initial forest ----------------------------------------------------------
initial_forest <- function(){
  pue_forest <- emptyforest(1)
  pue_forest$treeData$Species <- "Quercus ilex"
  pue_forest$treeData$DBH <- 9.11
  pue_forest$treeData$Height <- 530.2
  pue_forest$treeData$N <- 1750
  pue_forest$treeData$Z50 <- 200
  pue_forest$treeData$Z95 <- 1000
  pue_forest$treeData$LAI <- 2.2
  return(pue_forest)
}
forest <- initial_forest()

# Weather -----------------------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo
row.names(meteo) <- NULL
meteo <- rbind(meteo, meteo, meteo)
meteo$dates <- seq(as.Date("2001-01-01"), as.Date("2003-12-31"), by="day")
meteo$DOY <- 200
meteo$JulianDay <- 200
meteo$MinTemperature <- 20
meteo$MaxTemperature <- 30
meteo$MinRelativeHumidity <- 30
meteo$MaxRelativeHumidity <- 90
meteo$WindSpeed <- 2
meteo$Radiation  <- 30
meteo$Precipitation <- 0 # To force dessication


# Build initial input object ----------------------------------------------
if(model=="Sureau") {
  control <- defaultControl("Cochard")
  control$subdailyResults <- FALSE
  control$cavitationRefillStem <- "none"
  control$cavitationRefillLeaves <- "none"
  control$bareSoilEvaporation <- FALSE
  control$plantCapacitance <- TRUE
  control$cavitationFlux <- TRUE
  control$sapFluidityVariation <- TRUE
  control$leafCuticularTranspiration <- TRUE
  control$stemCuticularTranspiration <- TRUE
  control$rhizosphereOverlap <- "total"
  control$stomatalSubmodel <- "Baldocchi"
  control$sunlitShade <- FALSE
  control$gs_NightFrac <- 0.05
  control$verbose <- FALSE
} else {
  control <- defaultControl("Sperry")
  control$subdailyResults <- FALSE
  control$cavitationRefillStem <- "none"
  control$cavitationRefillLeaves <- "total"
  control$leafCavitationEffects <- FALSE
  control$bareSoilEvaporation <- FALSE
  control$sapFluidityVariation <- TRUE
  control$rhizosphereOverlap <- "total"
  control$sunlitShade <- FALSE
  control$verbose <- FALSE
}

x_initial <- forest2spwbInput(forest, soil, SpParamsMED, control)
x_initial$canopy$Tair <- 29
x_initial$canopy$Cair <- 386
x_initial$canopy$VPair <- 1.688
x_initial$soil$Temp <- c(32,29,27.71661)


# Test parameter modification ---------------------------------------------
tarcoh <- "T1_168"
if(model=="Sureau") {
  parNames = c("rfc@2",
               paste0(tarcoh,c("/LAI_live", "/Z95", "/Plant_kmax", "/Gswmax", "/Vmax298", "/Gs_P50", "/Gswmin", "/VC_P50", "/Vsapwood")))
  parMin <- c(10,1, 500, 0.3, 0.1, 20, -4.0, 0.001, -12.0, 1.0)
  parMax <- c(90,7, 6000, 4.0, 0.5, 100, -0.7, 0.010, -1.0, 30.0)
} else {
  parNames = c("rfc@2",
               paste0(tarcoh,c("/LAI_live", "/Z95", "/Plant_kmax", "/Gswmax", "/Vmax298", "/Gswmin", "/VC_d")))
  parMin <- c(10,1, 500, 0.3, 0.1, 20, 0.001, -12.0)
  parMax <- c(90,7, 6000, 4.0, 0.5, 100, 0.010, -1.0)
}
red_exclude <- c(1,2)
parMin_red <-parMin[-red_exclude]
parMax_red <-parMax[-red_exclude]
parNames_red <- parNames[-red_exclude]

customParams <- parMin
names(customParams) <- parNames
x_mod <- modifyInputParams(x_initial, customParams)

# Optimization functions --------------------------------------------------

# First day with maximum gs < 10% of gs initial in Sunlit leaves
sf_timetoclosure<-function(x){
  gswmax <- x$SunlitLeaves$GSWMax
  if(all(is.na(gswmax))) return(1)
  ttc <- which(gswmax < gswmax[1]*0.1)[1] 
  if(is.na(ttc)) return(length(gswmax))
  return(ttc)
}

#First day with PLC > 90% or missing StemPLC, meaning it stoped before
sf_timetofailure<-function(x){ 
  plc <- as.numeric(x$Plants$StemPLC[,1])
  # print(plc)
  na_plc <- which(is.na(plc))[1]
  n90 <- which(plc > 0.90)[1]
  n90 <- max(n90, sf_timetoclosure(x)) # Ensure time is not less than stomatal closure
  if(!is.na(na_plc)) return(min(n90, na_plc))
  return(n90)
}
# Survival time as the difference
sf_survivaltime<-function(x){ 
  max(0, sf_timetofailure(x) - sf_timetoclosure(x))
}

# Cumulated GPP to death
sf_gpp<-function(x){
  sum(x$Plants$GrossPhotosynthesis, na.rm=TRUE)
}
of_timetoclosure<-optimization_function(parNames = parNames,
                                     x = x_initial,
                                     meteo = meteo,
                                     latitude = pue_latitude, elevation = pue_elevation,
                                     slope = pue_slope, aspect = pue_aspect,
                                     summary_function = sf_timetoclosure)

of_timetoclosure_red<-optimization_function(parNames = parNames_red,
                                        x = x_initial,
                                        meteo = meteo,
                                        latitude = pue_latitude, elevation = pue_elevation,
                                        slope = pue_slope, aspect = pue_aspect,
                                        summary_function = sf_timetoclosure)
of_timetofailure<-optimization_function(parNames = parNames,
                                     x = x_initial,
                                     meteo = meteo,
                                     latitude = pue_latitude, elevation = pue_elevation,
                                     slope = pue_slope, aspect = pue_aspect,
                                     summary_function = sf_timetofailure)
of_timetofailure_red<-optimization_function(parNames = parNames_red,
                                        x = x_initial,
                                        meteo = meteo,
                                        latitude = pue_latitude, elevation = pue_elevation,
                                        slope = pue_slope, aspect = pue_aspect,
                                        summary_function = sf_timetofailure)

of_survivaltime<-optimization_function(parNames = parNames,
                                        x = x_initial,
                                        meteo = meteo,
                                        latitude = pue_latitude, elevation = pue_elevation,
                                        slope = pue_slope, aspect = pue_aspect,
                                        summary_function = sf_survivaltime)
of_survivaltime_red<-optimization_function(parNames = parNames_red,
                                       x = x_initial,
                                       meteo = meteo,
                                       latitude = pue_latitude, elevation = pue_elevation,
                                       slope = pue_slope, aspect = pue_aspect,
                                       summary_function = sf_survivaltime)

of_gpp<-optimization_function(parNames = parNames,
                                       x = x_initial,
                                       meteo = meteo,
                                       latitude = pue_latitude, elevation = pue_elevation,
                                       slope = pue_slope, aspect = pue_aspect,
                                       summary_function = sf_gpp)
of_gpp_red<-optimization_function(parNames = parNames_red,
                                           x = x_initial,
                                           meteo = meteo,
                                           latitude = pue_latitude, elevation = pue_elevation,
                                           slope = pue_slope, aspect = pue_aspect,
                                           summary_function = sf_gpp)
of_timetoclosure(parMin)
of_timetoclosure(parMax)
of_timetofailure(parMin)
of_timetofailure(parMax)
of_survivaltime(parMin)
of_survivaltime(parMax)
of_gpp(parMin)
of_gpp(parMax)


# Parameter matrices ------------------------------------------------------
p <- length(parMin)
M1 <- sweep(matrix(runif(p * n), nrow = n),2,STATS = (parMax-parMin), FUN = "*")
M1 <- sweep(M1, 2, parMin, FUN="+")
X1 <- data.frame(M1)
names(X1)<-parNames
M2 <- sweep(matrix(runif(p * n), nrow = n),2,STATS = (parMax-parMin), FUN = "*")
M2 <- sweep(M2, 2, parMin, FUN="+")
X2 <- data.frame(M2)
names(X2)<-parNames

X1_red <- X1[,-red_exclude]
X2_red <- X2[,-red_exclude]

# Parallelization functions --------------------------------------------------------
mult_timetoclosure <- function(X) {
  cat(paste0("Entering mult_timetoclosure n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_timetoclosure(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_timetofailure <- function(X) {
  cat(paste0("Entering mult_timetofailure n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_timetofailure(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_survivaltime <- function(X) {
  cat(paste0("Entering mult_survivaltime n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_survivaltime(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_gpp <- function(X) {
  cat(paste0("Entering mult_gpp n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_gpp(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_timetoclosure_red <- function(X) {
  cat(paste0("Entering mult_timetoclosure n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_timetoclosure_red(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_timetofailure_red <- function(X) {
  cat(paste0("Entering mult_timetofailure n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_timetofailure_red(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_survivaltime_red <- function(X) {
  cat(paste0("Entering mult_survivaltime n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_survivaltime_red(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}
mult_gpp_red <- function(X) {
  cat(paste0("Entering mult_gpp n = ", nrow(X), "\n"))
  doParallel::registerDoParallel(cores = ncores)
  r <- foreach::foreach(index = 1:nrow(X)) %dopar% {
    v <- as.numeric(X[index,])
    names(v) <- names(X)
    of_gpp_red(v, verbose = FALSE)
  }
  doParallel::stopImplicitCluster()
  r<- unlist(r)
  print(summary(r))
  return(r)
}

# mult_timetoclosure(X1[1:ncores,])
# mult_timetofailure(X1[1:ncores,])
# mult_survivaltime(X1[1:ncores,])
# mult_gpp(X1[1:ncores,])

# Sensitivity analyses ----------------------------------------------------

#FULL
sa_salt_timetoclosure <- sobolSalt(model = mult_timetoclosure, X1, X2,
                                   scheme="A", nboot = nboot)
saveRDS(sa_salt_timetoclosure, paste0("Rdata/sensitivity/",model, "/sa_salt_timetoclosure.rds"))

sa_salt_timetofailure <- sobolSalt(model = mult_timetofailure, X1, X2,
                                   scheme="A", nboot = nboot)
saveRDS(sa_salt_timetofailure, paste0("Rdata/sensitivity/", model, "/sa_salt_timetofailure.rds"))

sa_salt_survivaltime <- sobolSalt(model = mult_survivaltime, X1, X2,
                                   scheme="A", nboot = nboot)
saveRDS(sa_salt_survivaltime, paste0("Rdata/sensitivity/", model, "/sa_salt_survivaltime.rds"))

sa_salt_gpp <- sobolSalt(model = mult_gpp, X1, X2,
                                  scheme="A", nboot = nboot)
saveRDS(sa_salt_gpp, paste0("Rdata/sensitivity/", model, "/sa_salt_gpp.rds"))

#REDUCED
sa_salt_timetoclosure_red <- sobolSalt(model = mult_timetoclosure_red, X1_red, X2_red,
                                   scheme="A", nboot = nboot)
saveRDS(sa_salt_timetoclosure_red, paste0("Rdata/sensitivity/", model, "/sa_salt_timetoclosure_red.rds"))

sa_salt_timetofailure_red <- sobolSalt(model = mult_timetofailure_red, X1_red, X2_red,
                                   scheme="A", nboot = nboot)
saveRDS(sa_salt_timetofailure_red, paste0("Rdata/sensitivity/", model, "/sa_salt_timetofailure_red.rds"))

sa_salt_survivaltime_red <- sobolSalt(model = mult_survivaltime_red, X1_red, X2_red,
                                  scheme="A", nboot = nboot)
saveRDS(sa_salt_survivaltime_red, paste0("Rdata/sensitivity/", model, "/sa_salt_survivaltime_red.rds"))

sa_salt_gpp_red <- sobolSalt(model = mult_gpp_red, X1_red, X2_red,
                                      scheme="A", nboot = nboot)
saveRDS(sa_salt_gpp_red, paste0("Rdata/sensitivity/", model, "/sa_salt_gpp_red.rds"))
