# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_FontBlanche.R")


# Terrain -----------------------------------------------------------------
fb_latitude <- 43.24
fb_elevation <- 420
fb_slope <- 0 
fb_aspect <- 0

# Meteo -------------------------------------------------------------------
fb_meteo_raw <- read_delim("Data/FontBlanche/Climate_FontBlanche_GapFilled.csv",
                       delim = ";", escape_double = FALSE, trim_ws = TRUE)
fb_meteo <- read_delim("Data/FontBlanche/Climate_FontBlanche_GapFilled.csv",
                       delim = ";", escape_double = FALSE, trim_ws = TRUE) |>
  rename(dates = Date,
         MinTemperature = Tair_min,
         MaxTemperature = Tair_max,
         MeanTemperature = Tair_mean,
         Radiation = RG_sum,
         Precipitation = PPT_sum, 
         MinRelativeHumidity = RHair_min, 
         MaxRelativeHumidity = RHair_max, 
         MeanRelativeHumidity = RHair_mean,
         WindSpeed = WS_mean) |>
  mutate(dates = as.Date(dates, format = "%d/%m/%Y"),
         Radiation = Radiation/100)|>
  filter(dates >= as.Date("2017-01-01"),
         dates <= as.Date("2017-12-31"))



# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "annual"
control$cavitationRefillLeaves <- "annual"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE
control$sunlitShade <- FALSE

## Total overlap
control$rhizosphereOverlap <- "total"
x1t <- fontblanche_input(control)
S1t <- spwb(x1t, fb_meteo,
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1t, "Rdata/FontBlanche/Real_FontBlanche_Sperry_onepool.rds")


## Partial overlap
control$rhizosphereOverlap <- "partial"
x1p <- fontblanche_input(control)
S1p <- spwb(x1p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1p, "Rdata/FontBlanche/Real_FontBlanche_Sperry_partialpools.rds")



# Cochard initialization and run ----------------------------------------------------------
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
control$stomatalSubmodel <- "Baldocchi"

## Total overlap
control$rhizosphereOverlap <- "total"
x2t <- fontblanche_input(control)
S2t <- spwb(x2t, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2t, "Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_onepool.rds")


## Partial overlap
control$rhizosphereOverlap <- "partial"
x2p <- fontblanche_input(control)
S2p <- spwb(x2p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2p, "Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_partialpools.rds")


control$rhizosphereOverlap <- "total"
control$soilDisconnection <- TRUE
x2s <- fontblanche_input(control)
#Change canopy and soil variables
x2s$canopy$Tair <- 29
x2s$canopy$Cair <- 386
x2s$canopy$VPair <- 1.688
x2s$soil$Temp <- c(32,29,27.71661)

S2s <- spwb(x2s, fb_meteo,
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2s, "Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_disconnection.rds")


