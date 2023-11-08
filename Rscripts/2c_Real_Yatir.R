# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/Ancillary.R")

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
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- TRUE
control$rhizosphereOverlap <- "total"
control$sunlitShade <- FALSE

#Initialize input
x1 <- puechabon_input(control)

S1 <- spwb(x1, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S1, "Rdata/Yatir/Real_Yatir_Sperry.rds")

# shinyplot(S1)

# Cochard initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "annual"
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

S2 <- spwb(x2, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S2, "Rdata/Yatir/Real_Yatir_Sureau_Jarvis.rds")

control$stomatalSubmodel <- "Baldocchi"
control$leafCuticularTranspiration <- FALSE
x2b <- puechabon_input(control)

#Call simulation function
S2b <- spwb(x2b, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S2b, "Rdata/Yatir/Real_Yatir_Sureau_Baldocchi.rds")

# Sperry (segmented) ------------------------------------------------------
x3 <- x1
wb <- hydraulics_psi2Weibull(psi50 = -4.0, psi88 = -4.5)
x3$paramsTranspiration$VCleaf_c <- wb["c"]
x3$paramsTranspiration$VCleaf_d <- wb["d"]
x3$control$leafCavitationEffects <- FALSE
S3 <- spwb(x3, yat_meteo, 
           latitude = yat_latitude, elevation = yat_elevation, 
           slope = yat_slope, aspect = yat_aspect)
saveRDS(S3, "Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")

# Plots -------------------------------------------------------------------
# p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
# p2 <- plot(S2, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
# p3 <- plot(S3, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-seg")+theme(legend.position = c(0.8,0.8))
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/SoilPsi_Yatir_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
# p2 <- plot(S2, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
# p3 <- plot(S3, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/HydraulicRedistribution_Yatir_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/LeafPsiRange_Yatir_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/Transpiration_Yatir_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/StomatalConductance_Yatir_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/StemPLC_Yatir_Real.png", p, width = 8, height = 11)
# 
# 
# p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/LeafPLC_Yatir_Real.png", p, width = 8, height = 11)
# 
# p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
# p2 <- plot(S2, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
# p3 <- plot(S3, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-seg")
# p <-plot_grid(p1, p3, p2,
#               nrow = 3)
# ggsave2("Plots/Yatir_Real/SoilPlantConductance_Yatir_Real.png", p, width = 8, height = 11)

# shinyplot(S2)
