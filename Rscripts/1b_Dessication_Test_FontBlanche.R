# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_FontBlanche.R")

# Terrain -----------------------------------------------------------------
fb_latitude <- 43.24
fb_elevation <- 420
fb_slope <- 0
fb_aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo[1:180,]
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
control$cavitationRefillStem <- "none"
control$cavitationRefillLeaves <- "none"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE
control$sunlitShade <- FALSE

#Initialize input
control$rhizosphereOverlap <- "total"
x1t <- fontblanche_input(control)
#Change canopy and soil variables
x1t$canopy$Tair <- 29
x1t$canopy$Cair <- 386
x1t$canopy$VPair <- 1.688
x1t$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1t <- spwb(x1t, meteo, 
           latitude = fb_latitude, elevation = fb_elevation, 
           slope = fb_slope, aspect = fb_aspect)
saveRDS(S1t, "Rdata/FontBlanche/Dessication_FontBlanche_Sperry_onepool.rds")

#Initialize input
control$rhizosphereOverlap <- "partial"
x1p <- fontblanche_input(control)
#Change canopy and soil variables
x1p$canopy$Tair <- 29
x1p$canopy$Cair <- 386
x1p$canopy$VPair <- 1.688
x1p$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1p <- spwb(x1p, meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1p, "Rdata/FontBlanche/Dessication_FontBlanche_Sperry_partialpools.rds")



# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "none"
control$cavitationRefillLeaves <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$cavitationFlux <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- FALSE
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001

control$rhizosphereOverlap <- "total"
x2t <- fontblanche_input(control)
#Change canopy and soil variables
x2t$canopy$Tair <- 29
x2t$canopy$Cair <- 386
x2t$canopy$VPair <- 1.688
x2t$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2t <- spwb(x2t, meteo, 
           latitude = fb_latitude, elevation = fb_elevation, 
           slope = fb_slope, aspect = fb_aspect)
saveRDS(S2t, "Rdata/FontBlanche/Dessication_FontBlanche_Sureau_Jarvis_onepool.rds")


control$stomatalSubmodel <- "Baldocchi"
x2b <- fontblanche_input(control)

#Change canopy and soil variables
x2b$canopy$Tair <- 29
x2b$canopy$Cair <- 386
x2b$canopy$VPair <- 1.688
x2b$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2b <- spwb(x2b, meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2b, "Rdata/FontBlanche/Dessication_FontBlanche_Sureau_Baldocchi_onepool.rds")

control$rhizosphereOverlap <- "partial"
x2p <- fontblanche_input(control)
#Change canopy and soil variables
x2p$canopy$Tair <- 29
x2p$canopy$Cair <- 386
x2p$canopy$VPair <- 1.688
x2p$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2p <- spwb(x2p, meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2p, "Rdata/FontBlanche/Dessication_FontBlanche_Sureau_Baldocchi_partialpools.rds")


control$rhizosphereOverlap <- "total"
control$soilDisconnection <- TRUE
x2s <- fontblanche_input(control)
#Change canopy and soil variables
x2s$canopy$Tair <- 29
x2s$canopy$Cair <- 386
x2s$canopy$VPair <- 1.688
x2s$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2s <- spwb(x2s, meteo,
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2s, "Rdata/FontBlanche/Dessication_FontBlanche_Sureau_Baldocchi_disconnection.rds")

