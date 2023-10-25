# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/Ancillary.R")

FrFBn_Qilex_LAI2 <- read_table("Data/Puechabon/FrFBn_Qilex_LAI2.csv")


# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270

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
control$subdailyResults <- FALSE
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE

#Initialize input
x1 <- puechabon_input(control)

S1 <- spwb(x1, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
# shinyplot(S1)

# Cochard initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- FALSE
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE
control$sapFluidityVariation <- FALSE
control$rhizosphereOverlap <- "partial"

#Initialize input
x2 <- puechabon_input(control)

S2 <- spwb(x2, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)

# Plots -------------------------------------------------------------------
p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")
p2 <- plot(S2, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Real/SoilPsi_Puechabon_Real.png", p, width = 6, height = 8)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Real/LeafPsiRange_Puechabon_Real.png", p, width = 6, height = 8)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,2))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "Transpiration", bySpecies = TRUE)+ylim(c(0,2))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Real/Transpiration_Puechabon_Real.png", p, width = 6, height = 8)


p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Real/StomatalConductance_Puechabon_Real.png", p, width = 6, height = 8)

# shinyplot(S2)