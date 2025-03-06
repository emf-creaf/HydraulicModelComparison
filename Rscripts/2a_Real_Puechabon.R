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
control$stemCavitationRecovery <- "annual"
control$leafCavitationRecovery <- "annual"
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
control$stomatalSubmodel <- "Jarvis"
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001

#Initialize input
control$stomatalSubmodel <- "Baldocchi"
control$leafCuticularTranspiration <- FALSE
x2b <- puechabon_input(control)

#Call simulation function
S2b <- spwb(x2b, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
saveRDS(S2b, "Rdata/Puechabon/Real_Puechabon_Sureau_Baldocchi.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
# gs_psi50 <- x2b$paramsTranspiration$Gs_P50
# gs_slope <- x2b$paramsTranspiration$Gs_slope
# gs_psi88 <- gs_psi50 + log((1/0.88)-1.0)*(25.0/gs_slope)
# wb <- hydraulics_psi2Weibull(psi50 = gs_psi50, psi88 = gs_psi88)
# From Limousin et al. (2022)
# Quercus ilex Kleaf P50 = -5 P88 = -6.5
wb <- hydraulics_psi2Weibull(psi50 = -5.0, psi88 = -6.5)
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

medfate:::.updateBelow(x1s)

S1s <- spwb(x1s, pue_meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)
saveRDS(S1s, "Rdata/Puechabon/Real_Puechabon_Sperry_segmented.rds")

# Plots -------------------------------------------------------------------
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry.rds")
S2b <- readRDS("Rdata/Puechabon/Real_Puechabon_Sureau_Baldocchi.rds")
S1s <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry_segmented.rds")
p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/SoilPsi_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry")
p2 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p3 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/HydraulicRedistribution_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/LeafPsiRange_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/Transpiration_Puechabon_Real.png", p, width = 8, height = 11)


p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/StomatalConductance_Puechabon_Real.png", p, width = 8, height = 11)


p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/StemPLC_Puechabon_Real.png", p, width = 8, height = 11)


p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/LeafPLC_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p3 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/SoilPlantConductance_Puechabon_Real.png", p, width = 8, height = 11)


library(readxl)
library(readr)
E_data <- read_excel("Data/Puechabon/Sapflow_day_SMR_2016-2018.xlsx") |>
  dplyr::select(Date, E_Av.gap)|>
  dplyr::rename(dates = Date, 
                E_T1_168 = E_Av.gap) |>
  dplyr::mutate(dates = as.Date(dates),
                E_T1_168 = E_T1_168/x1$above$LAI_live) |>
  as.data.frame()
row.names(E_data) <- as.character(E_data$dates)

wp_data <- read.table("Data/Puechabon/Water_Potential_MIND_Control-1.csv",sep = ";", 
                      dec=",", header = TRUE) |>
  dplyr::mutate(Date = as.Date(Date, "%d/%m/%Y")) |>
  dplyr::filter(Treatment == "Control") |>
  dplyr::select(-Treatment)  |>
  dplyr::group_by(Date) |>
  summarise(PD_T1_168 = mean(Pd_pot, na.rm=TRUE),
            PD_T1_168_err = sqrt(var(Pd_pot, na.rm=TRUE)),
            MD_T1_168 = mean(Md_pot, na.rm=TRUE),
            MD_T1_168_err = sqrt(var(Md_pot, na.rm=TRUE))) |>
  as.data.frame()
row.names(wp_data) <- as.character(wp_data$Date)
wp_evaluation <- function(S, wp_data, E_data, title) {
  p1 <- evaluation_plot(S, E_data, type="E", cohort = "T1_168")+ ylim(c(0,2))+
    labs(title = title, subtitle = "Sap flux Quercus ilex")+
    theme_classic()+theme(legend.position = c(0.05,0.8))
  p2<- evaluation_plot(S, wp_data, type="WP", cohort = "T1_168")+ ylim(c(-10,0))+
    labs(subtitle = "Leaf WP Quercus ilex", title = "")+theme_classic()+ theme(legend.position = "none")
  p3 <- plot(S, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+
    theme_classic()+ theme(legend.position = "none")+labs(title=title, subtitle="PLC")
  return(cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(2,2,1.3)))  
}
p1 <-wp_evaluation(S1, wp_data, E_data, "Sperry non-segmented")
p2 <-wp_evaluation(S1s, wp_data, E_data,"Sperry segmented")
p3 <-wp_evaluation(S2b, wp_data, E_data,"Sureau")
p <-plot_grid(p1, p2, p3, nrow = 3)
ggsave2("Plots/Puechabon_Real/Evaluation_Puechabon_Real.png", p, width = 20, height = 10, bg = "white")
