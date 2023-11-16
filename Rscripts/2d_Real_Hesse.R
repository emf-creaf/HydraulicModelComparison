# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0d_Ancillary_Hesse.R")

# Terrain -----------------------------------------------------------------
hes_latitude <- 48.6
hes_elevation <- 300
hes_slope <- 0
hes_aspect <- 0

# Soil --------------------------------------------------------------------
hes_soil <- hesse_soil()
sum(soil_waterFC(hes_soil, "VG"))
sum(soil_waterExtractable(hes_soil, "VG", -6))


# Meteo -------------------------------------------------------------------
# hes_meteo <- read.table("Data/Hesse/FRAHES_meteoData.txt", 
#                         sep = "\t", header = TRUE, dec=".") |>
#   rename(dates = Date) |>
#   mutate(dates = as.Date(dates))

hes_meteo <- read.table("Data/Hesse/clim_Hesse.txt", 
                        sep = "\t", header = TRUE, dec=".") |>
  mutate(prec = prec*10) |>
  rename(MinTemperature = tmin,
         MaxTemperature = tmax,
         MeanTemperature = tmean,
         Radiation = glo,
         Precipitation = prec, 
         MinRelativeHumidity = RHmin, 
         MaxRelativeHumidity = RHmax,
         WindSpeed = wnd) |>
  mutate(dates = seq(as.Date("1999-01-01"), as.Date("2010-12-31"), by="day")) |>
  filter(dates >= as.Date("2009-01-01"),
         dates <= as.Date("2010-12-31"))


# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "annual"
control$cavitationRefillLeaves <- "annual"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- TRUE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE

#Initialize input
x1 <- hesse_input(control)

S1 <- spwb(x1, hes_meteo, latitude = hes_latitude, elevation = hes_elevation)
saveRDS(S1, "Rdata/Hesse/Real_Hesse_Sperry.rds")

# shinyplot(S1)

# Cochard initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "annual"
control$cavitationRefillLeaves <- "annual"
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

#Initialize input
x2 <- hesse_input(control)

S2 <- spwb(x2, hes_meteo, latitude = hes_latitude, elevation = hes_elevation)
saveRDS(S2, "Rdata/Hesse/Real_Hesse_Sureau_Jarvis.rds")

control$stomatalSubmodel <- "Baldocchi"
control$leafCuticularTranspiration <- FALSE
x2b <- hesse_input(control)

#Call simulation function
S2b <- spwb(x2b, hes_meteo, latitude = hes_latitude, elevation = hes_elevation)
saveRDS(S2b, "Rdata/Hesse/Real_Hesse_Sureau_Baldocchi.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$cavitationRefillLeaves <- "total"
x1s$control$leafCavitationEffects <- FALSE
gs_psi50 <- x2j$paramsTranspiration$Gs_P50
gs_slope <- x2j$paramsTranspiration$Gs_slope
gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
S1s <- spwb(x1s, hes_meteo, 
           latitude = hes_latitude, elevation = hes_elevation, 
           slope = hes_slope, aspect = hes_aspect)
saveRDS(S1s, "Rdata/Hesse/Real_Hesse_Sperry_segmented.rds")
