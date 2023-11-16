# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0a_Ancillary_Puechabon.R")

FrFBn_Qilex_LAI2 <- read_table("Data/Puechabon/FrFBn_Qilex_LAI2.csv")


# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270
pue_slope <- 0
pue_aspect <- 0

# Soil --------------------------------------------------------------------
pue_soil <- puechabon_soil()
sum(soil_waterFC(pue_soil, "VG"))
sum(soil_waterExtractable(pue_soil, "VG", -6))


# Meteo -------------------------------------------------------------------
pue_meteo <- read_delim("Data/Puechabon/Climat_Puechabon_site.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE) |>
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
  filter(dates >= as.Date("2016-01-01"),
         dates <= as.Date("2018-12-31"))



# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "annual"
control$cavitationRefillLeaves <- "annual"
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- TRUE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE

#Initialize input
x1 <- puechabon_input(control)

S1 <- spwb(x1, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
saveRDS(S1, "Rdata/Puechabon/Real_Puechabon_Sperry.rds")

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
x2 <- puechabon_input(control)

S2 <- spwb(x2, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
saveRDS(S2, "Rdata/Puechabon/Real_Puechabon_Sureau_Jarvis.rds")

control$stomatalSubmodel <- "Baldocchi"
control$leafCuticularTranspiration <- FALSE
x2b <- puechabon_input(control)

#Call simulation function
S2b <- spwb(x2b, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
saveRDS(S2b, "Rdata/Puechabon/Real_Puechabon_Sureau_Baldocchi.rds")

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
medfate:::.updateBelow(x1s)

S1s <- spwb(x1s, pue_meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)
saveRDS(S1s, "Rdata/Puechabon/Real_Puechabon_Sperry_segmented.rds")

# Plots -------------------------------------------------------------------
# p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
# p2 <- plot(S2, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
# p3 <- plot(S3, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-seg")+theme(legend.position = c(0.8,0.8))
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/SoilPsi_Puechabon_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
# p2 <- plot(S2, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
# p3 <- plot(S3, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/HydraulicRedistribution_Puechabon_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/LeafPsiRange_Puechabon_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/Transpiration_Puechabon_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/StomatalConductance_Puechabon_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/StemPLC_Puechabon_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/LeafPLC_Puechabon_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Puechabon_Real/SoilPlantConductance_Puechabon_Real.png", p, width = 8, height = 11)

# shinyplot(S2)
