# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_FontBlanche.R")
source("Rscripts/0_ExtractParams.R")

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
control$stemCavitationRecovery <- "none"
control$leafCavitationRecovery <- "none"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE
control$sunlitShade <- FALSE

#Initialize input
control$rhizosphereOverlap <- "total"
x1 <- fontblanche_input(control)
#Change canopy and soil variables
x1$canopy$Tair <- 29
x1$canopy$Cair <- 386
x1$canopy$VPair <- 1.688
x1$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1 <- spwb(x1, meteo, 
           latitude = fb_latitude, elevation = fb_elevation, 
           slope = fb_slope, aspect = fb_aspect)
saveRDS(S1, "Rdata/FontBlanche/Dessication_FontBlanche_Sperry.rds")

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
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001

control$rhizosphereOverlap <- "total"
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


# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$rhizosphereOverlap <- "total"
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
# gs_psi50 <- x2b$paramsTranspiration$Gs_P50
# gs_slope <- x2b$paramsTranspiration$Gs_slope
# gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
# From Limousin et al. (2022)
# Quercus ilex Kleaf P50 = -5 P88 = -6.5
# P. halepensis P50 = -2 P88 = -2.2
wb_1 <- hydraulics_psi2Weibull(psi50 = -2.0, psi88 = -2.2)
wb_2 <- hydraulics_psi2Weibull(psi50 = -5.0, psi88 = -6.5)
x1s$paramsTranspiration$VCleaf_c <- c(wb_1["c"], wb_2["c"]) 
x1s$paramsTranspiration$VCleaf_d <- c(wb_1["d"], wb_2["d"])
# Bartlett et al. (2016)
# Root P50 = 0.4892 + 0.742 * Stem_P50
root_P50 <- 0.4892 + 0.742*x2b$paramsTranspiration$VCstem_P50
# P88,stem=−1.4264+1.2593* P50,stem
root_P88 <- -1.4264 + 1.2593*root_P50
wb_1 <- hydraulics_psi2Weibull(psi50 = root_P50[1], psi88 = root_P88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = root_P50[2], psi88 = root_P88[2])
x1s$paramsTranspiration$VCroot_c <- c(wb_1["c"], wb_2["c"]) 
x1s$paramsTranspiration$VCroot_d <- c(wb_1["d"], wb_2["d"])

medfate:::.updateBelow(x1s)
S1s <- spwb(x1s, meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1s, "Rdata/FontBlanche/Dessication_FontBlanche_Sperry_segmented.rds")


# Plots -------------------------------------------------------------------
library(medfate)
library(cowplot)
library(ggplot2)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/FontBlanche/Dessication_FontBlanche_Sperry.rds")
S1s <- readRDS("Rdata/FontBlanche/Dessication_FontBlanche_Sperry_segmented.rds")
S2b <- readRDS("Rdata/FontBlanche/Dessication_FontBlanche_Sureau_Baldocchi_disconnection.rds")


p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/SoilPsi_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
p2 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p3 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/HydraulicRedistribution_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/LeafPsiRange_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/Transpiration_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/StomatalConductance_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/StemPLC_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/LeafPLC_FontBlanche_Dessication.png", p, width = 6, height = 11)

p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/FontBlanche_Dessication/SoilPlantConductance_FontBlanche_Dessication.png", p, width = 6, height = 11)

extract_spparams(x2b)
