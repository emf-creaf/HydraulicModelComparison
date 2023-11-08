# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(meteoland)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/0_Ancillary.R")

# Terrain -----------------------------------------------------------------
hes_latitude <- 48.6
hes_elevation <- 300
hes_slope <- 0
hes_aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo[1:50,]
meteo$DOY <- 200
meteo$JulianDay <- 200
meteo$MinTemperature <- 20
meteo$MaxTemperature <- 30
meteo$MinRelativeHumidity <- 30
meteo$MaxRelativeHumidity <- 90
meteo$WindSpeed <- 2
meteo$Radiation  <- 30
meteo$Precipitation <- 0 # To force dessication

# Sperry simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- TRUE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE

#Initialize input
x1 <- hesse_input(control)

#Change canopy and soil variables
x1$canopy$Tair <- 29
x1$canopy$Cair <- 386
x1$canopy$VPair <- 1.688
x1$soil$Temp <- c(32,29,27.71661)
x1$paramsInterception

#Call simulation function
x1$internalPhenology$phi <- 1 # Unfolded leaves
x1$internalPhenology$gdd <- 2000

S1 <- spwb(x1, meteo, 
           latitude = hes_latitude, elevation = hes_elevation, 
           slope = hes_slope, aspect = hes_aspect)

saveRDS(S1, "Rdata/Hesse/Dessication_Hesse_Sperry.rds")

# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$cavitationFlux <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- FALSE
control$rhizosphereOverlap <- "total"
control$stomatalSubmodel <- "Jarvis"
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001

x2j <- hesse_input(control)
x2j$internalPhenology$phi <- 1 # Unfolded leaves
x2j$internalPhenology$gdd <- 2000

#Change canopy and soil variables
x2j$canopy$Tair <- 29
x2j$canopy$Cair <- 386
x2j$canopy$VPair <- 1.688
x2j$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2j <- spwb(x2j, meteo, 
           latitude = hes_latitude, elevation = hes_elevation, 
           slope = hes_slope, aspect = hes_aspect)
saveRDS(S2j, "Rdata/Hesse/Dessication_Hesse_Sureau_Jarvis.rds")


control$stomatalSubmodel <- "Baldocchi"
x2b <- hesse_input(control)
x2b$internalPhenology$phi <- 1 # Unfolded leaves
x2b$internalPhenology$gdd <- 2000

#Change canopy and soil variables
x2b$canopy$Tair <- 29
x2b$canopy$Cair <- 386
x2b$canopy$VPair <- 1.688
x2b$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2b <- spwb(x2b, meteo, 
           latitude = hes_latitude, elevation = hes_elevation, 
           slope = hes_slope, aspect = hes_aspect)
saveRDS(S2b, "Rdata/Hesse/Dessication_Hesse_Sureau_Baldocchi.rds")


# Sperry (segmented) ------------------------------------------------------
x1s <- x1
gs_psi50 <- x2j$paramsTranspiration$Gs_P50
gs_slope <- x2j$paramsTranspiration$Gs_slope
gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)

wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
x1s$control$leafCavitationEffects <- FALSE
S1s <- spwb(x1s, meteo, 
           latitude = hes_latitude, elevation = hes_elevation, 
           slope = hes_slope, aspect = hes_aspect)
saveRDS(S1s, "Rdata/Hesse/Dessication_Hesse_Sperry_segmented.rds")


