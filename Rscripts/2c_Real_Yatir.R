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
control$stemCavitationRecovery <- "annual"
control$leafCavitationRecovery <- "annual"
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

# Sureau initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sureau")
control$subdailyResults <- TRUE
control$stemCavitationRecovery <- "annual"
control$leafCavitationRecovery <- "annual"
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
control$leafCuticularTranspiration <- FALSE
x2b <- yatir_input(control)

#Call simulation function
S2b <- spwb(x2b, yat_meteo, latitude = yat_latitude, elevation = yat_elevation)
saveRDS(S2b, "Rdata/Yatir/Real_Yatir_Sureau_Baldocchi.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
# gs_psi50 <- x2b$paramsTranspiration$Gs_P50
# gs_slope <- x2b$paramsTranspiration$Gs_slope
# gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
# wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
wb <- hydraulics_psi2Weibull(psi50 = -2.0, psi88 = -2.2)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
# Bartlett et al. (2016)
# Root P50 = 0.4892 + 0.742 * Stem_P50
root_P50 <- 0.4892 + 0.742*x2b$paramsTranspiration$VCstem_P50
# P88,stem=−1.4264+1.2593* P50,stem
root_P88 <- -1.4264 + 1.2593*root_P50
wb <- hydraulics_psi2Weibull(psi50 = root_P50, psi88 = root_P88)
x1s$paramsTranspiration$VCroot_c <- wb["c"] 
x1s$paramsTranspiration$VCroot_d <- wb["d"]

S1s <- spwb(x1s, yat_meteo, 
           latitude = yat_latitude, elevation = yat_elevation, 
           slope = yat_slope, aspect = yat_aspect)
saveRDS(S1s, "Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")

# Sperry segmented (TLP) --------------------------------------------------
x1st <- x1s
wb <- hydraulics_psi2Weibull(psi50 = -2.7, psi88 = -1.4264 + 1.2593*-2.7)
x1st$paramsTranspiration$VCleaf_c <- wb["c"]
x1st$paramsTranspiration$VCleaf_d <- wb["d"]
medfate:::.updateBelow(x1st)

S1st <- spwb(x1st, yat_meteo, 
            latitude = yat_latitude, elevation = yat_elevation, 
            slope = yat_slope, aspect = yat_aspect)
saveRDS(S1st, "Rdata/Yatir/Real_Yatir_Sperry_segmented_TLP.rds")

extract_spparams(x1)
extract_spparams(x1s)
extract_spparams(x2b)

# Plots -------------------------------------------------------------------
library(ggplot2)
library(cowplot)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/Yatir/Real_Yatir_Sperry.rds")
S1s <- readRDS("Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")
S1st <- readRDS("Rdata/Yatir/Real_Yatir_Sperry_segmented_TLP.rds")
S2b <- readRDS("Rdata/Yatir/Real_Yatir_Sureau_Baldocchi.rds")

p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented-tlp")+theme(legend.position = c(0.8,0.8))
p4 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/SoilPsi_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
p2 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p3 <- plot(S1st, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/HydraulicRedistribution_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/LeafPsiRange_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/Transpiration_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/StomatalConductance_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/StemPLC_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/LeafPLC_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S1st, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented-tlp")
p4 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/SoilPlantConductance_Yatir_Real.png", p, width = 8, height = 11)


wp_evaluation <- function(S,title) {
  p1 <- plot(S, "Transpiration", bySpecies = TRUE)+ylim(c(0,1.0))+
    labs(title = title, subtitle = "Transpiration Pinus halepensis")+
    theme_classic()+theme(legend.position = "none")
  p2<- plot(S, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-6.5,0))+
    labs(subtitle = "Leaf WP Pinus halepensis", title = title)+
    theme_classic()+ theme(legend.position = "none")
  p3 <- plot(S, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+
    theme_classic()+ theme(legend.position = "none")+
    labs(title=title, subtitle="PLC")
  return(cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(2,2,1.3)))  
}
p1 <-wp_evaluation(S1, "Sperry non-segmented")
p2 <-wp_evaluation(S1s,"Sperry segmented (Kleaf)")
p3 <-wp_evaluation(S1st,"Sperry segmented (TLP)")
p4 <-wp_evaluation(S2b,"Sureau")
p <-plot_grid(p1, p2, p3, p4, nrow = 4)
ggsave2("Plots/Yatir_Real/Evaluation_Yatir_Real.png", p, width = 20, height = 12, bg = "white")

