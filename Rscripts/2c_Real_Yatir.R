# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0c_Ancillary_Yatir.R")

# Terrain -----------------------------------------------------------------
yat_latitude <- 31.3
yat_elevation <- 650
yat_slope <- 0
yat_aspect <- 0

# Soil --------------------------------------------------------------------
yat_soil <- yatir_soil()
sum(soil_waterFC(yat_soil, "VG"))
sum(soil_waterExtractable(yat_soil, "VG", -6))


# Meteo -------------------------------------------------------------------
yat_meteo <- read.table("Data/Yatir/Yatir_Climate_ForSurEau_2010_2022.csv", 
                        sep = ";", header = TRUE, dec=",") |>
  rename(dates = DATE,
         MinTemperature = Tair_min,
         MaxTemperature = Tair_max,
         MeanTemperature = Tair_mean,
         Radiation = RG_sum,
         Precipitation = PPT_sum, 
         MinRelativeHumidity = RHair_min, 
         MaxRelativeHumidity = RHair_max, 
         MeanRelativeHumidity = RHair_mean,
         WindSpeed = WS_mean) |>
  mutate(dates = as.Date(dates, format = "%d/%m/%Y")) |>
  filter(dates >= as.Date("2020-01-01"),
         dates <= as.Date("2022-12-31"))



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
x1 <- yatir_input(control)

S1 <- spwb(x1, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S1, "Rdata/Yatir/Real_Yatir_Sperry.rds")

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
x2j <- yatir_input(control)

S2j <- spwb(x2j, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S2j, "Rdata/Yatir/Real_Yatir_Sureau_Jarvis.rds")

control$stomatalSubmodel <- "Baldocchi"
control$leafCuticularTranspiration <- FALSE
x2b <- yatir_input(control)

#Call simulation function
S2b <- spwb(x2b, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S2b, "Rdata/Yatir/Real_Yatir_Sureau_Baldocchi.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$cavitationRefillLeaves <- "total"
x1s$control$leafCavitationEffects <- FALSE
gs_psi50 <- x2b$paramsTranspiration$Gs_P50
gs_slope <- x2b$paramsTranspiration$Gs_slope
gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
S1s <- spwb(x1s, yat_meteo, 
           latitude = yat_latitude, elevation = yat_elevation, 
           slope = yat_slope, aspect = yat_aspect)
saveRDS(S1s, "Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")
