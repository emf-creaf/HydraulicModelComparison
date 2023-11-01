# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/Ancillary.R")

# k_RSApo <- k_RSApo1 + k_RSApo2 + k_RSApo3
# k_RSApo <- 2
# k_SLApo <- 4
# k_LSym <- 2
# 1/(1/k_RSApo + 1/k_SLApo + 1/k_LSym)


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



# Granier initialization and run ------------------------------------------
control <- defaultControl("Granier")
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$rhizosphereOverlap <- "total"

x0t <- fontblanche_input(control)
S0t <- spwb(x0t, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)

control$rhizosphereOverlap <- "partial"
x0p <- fontblanche_input(control)
#Call simulation function
S0p <- spwb(x0p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)

# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- FALSE
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE

## Total overlap
control$rhizosphereOverlap <- "total"
x1t <- fontblanche_input(control)
S1t <- spwb(x1t, fb_meteo,
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)

## Partial overlap
control$rhizosphereOverlap <- "partial"
x1p <- fontblanche_input(control)
S1p <- spwb(x1p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)



# Cochard initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- FALSE
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE
control$sapFluidityVariation <- FALSE

## Total overlap
control$rhizosphereOverlap <- "total"
x2t <- fontblanche_input(control)
S2t <- spwb(x2t, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)

## Partial overlap
control$rhizosphereOverlap <- "partial"
x2p <- fontblanche_input(control)
S2p <- spwb(x2p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)




# Plots -------------------------------------------------------------------
p1 <- plot(S1t, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (total)")
p2 <- plot(S1p, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (total)")
p4 <- plot(S2p, "SoilPsi")+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/SoilPsi_FontBlanche_Real.png", p, width = 10, height = 8)

p1 <- plot(S1t, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "HydraulicRedistribution")+ylim(c(0,1))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/HydraulicRedistribution_FontBlanche_Real.png", p, width = 10, height = 8)

p1 <- plot(S1t, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (total)")
p2 <- plot(S1p, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (total)")
p4 <- plot(S2p, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = c(0.2,0.2))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/LeafPsiRange_FontBlanche_Real.png", p, width = 10, height = 8)


p1 <- plot(S1t, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/StemPLC_FontBlanche_Real.png", p, width = 10, height = 8)


p1 <- plot(S1t, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.8,0.2))+labs(title="Sperry (total)")
p2 <- plot(S1p, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.8,0.2))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.2,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/LeafPLC_FontBlanche_Real.png", p, width = 10, height = 8)

p1 <- plot(S1t, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "Transpiration", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/Transpiration_FontBlanche_Real.png", p, width = 10, height = 8)

p1 <- plot(S1t, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (total)")
p2 <- plot(S1p, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry (partial)")
p3 <- plot(S2t, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (total)")
p4 <- plot(S2p, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,3))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau (partial)")
p <-plot_grid(p1, p2, p3, p4,
              nrow = 2)
ggsave2("Plots/FontBlanche_Real/SoilPlantConductance_FontBlanche_Real.png", p, width = 10, height = 8)

