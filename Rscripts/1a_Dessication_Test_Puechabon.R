# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(meteoland)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/0a_Ancillary_Puechabon.R")

# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270
pue_slope <- 0
pue_aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo[1:150,]
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
control$stemCavitationRecovery <- "none"
control$leafCavitationRecovery <- "none"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE

#Initialize input
x1 <- puechabon_input(control)

#Change canopy and soil variables
x1$canopy$Tair <- 29
x1$canopy$Cair <- 386
x1$canopy$VPair <- 1.688
x1$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1 <- spwb(x1, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)

saveRDS(S1, "Rdata/Puechabon/Dessication_Puechabon_Sperry.rds")

# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sureau")
control$subdailyResults <- TRUE
control$stemCavitationRecovery <- "none"
control$leafCavitationRecovery <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$cavitationFlux <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- FALSE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001
control$stomatalSubmodel <- "Baldocchi"
x2b <- puechabon_input(control)

#Change canopy and soil variables
x2b$canopy$Tair <- 29
x2b$canopy$Cair <- 386
x2b$canopy$VPair <- 1.688
x2b$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2b <- spwb(x2b, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)
saveRDS(S2b, "Rdata/Puechabon/Dessication_Puechabon_Sureau_Baldocchi.rds")


# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
gs_psi50 <- x2b$paramsTranspiration$Gs_P50
gs_slope <- x2b$paramsTranspiration$Gs_slope
gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
medfate:::.updateBelow(x1s)
S1s <- spwb(x1s, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)
saveRDS(S1s, "Rdata/Puechabon/Dessication_Puechabon_Sperry_segmented.rds")


# Plots -------------------------------------------------------------------
Sys.setlocale("LC_ALL", "en_US.utf8")

p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/SoilPsi_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
p2 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p3 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/HydraulicRedistribution_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/LeafPsiRange_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/Transpiration_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/StomatalConductance_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/StemPLC_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/LeafPLC_Puechabon_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Dessication/SoilPlantConductance_Puechabon_Dessication.png", p, width = 6, height = 11)

