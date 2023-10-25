# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/Ancillary.R")

# Terrain -----------------------------------------------------------------
fb_latitude <- 43.24
fb_elevation <- 420
fb_slope <- 0
fb_aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo[1:300,]
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



# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$sapFluidityVariation <- FALSE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE

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


# Plots -------------------------------------------------------------------
p1 <- plot(S1t, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (total)")
p2 <- plot(S1p, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (total)")
p4 <- plot(S2p, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Dessication/SoilPsi_FontBlanche_Dessication.png", p, width = 10, height = 8)

p1 <- plot(S1t, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
          nrow = 2)
ggsave2("Plots/FontBlanche_Dessication/LeafPsiRange_FontBlanche_Dessication.png", p, width = 10, height = 8)


p1 <- plot(S1t, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Dessication/StemPLC_FontBlanche_Dessication.png", p, width = 10, height = 8)


p1 <- plot(S1t, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.8,0.2))+labs(title="Sperry (total)")
p2 <- plot(S1p, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.8,0.2))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Dessication/LeafPLC_FontBlanche_Dessication.png", p, width = 10, height = 8)

p1 <- plot(S1t, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Dessication/Transpiration_FontBlanche_Dessication.png", p, width = 10, height = 8)

