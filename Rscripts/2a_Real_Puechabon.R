# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0a_Ancillary_Puechabon.R")
source("Rscripts/0_ExtractParams.R")

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
control$leafCavitationEffects <- FALSE
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
# From Limousin et al. (2022)
# Quercus ilex Kleaf P50 = -4.5 P88 = -7.0
psi50 = -4.5
slope = 20
psi88 <- psi50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb <- hydraulics_psi2Weibull(psi50 = psi50, psi88 = psi88)
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
psi12 <- medfate::hydraulics_xylemPsi(0.88,1, wb["c"], wb["d"])
x1s$paramsTranspiration$VCleaf_P50 <- psi50
x1s$paramsTranspiration$VCleaf_slope <- slope
# Bartlett et al. (2016)
# Root P50 = 0.4892 + 0.742 * Stem_P50
root_P50 <- 0.4892 + 0.742*x2b$paramsTranspiration$VCstem_P50
# P88,stem=−1.4264+1.2593* P50,stem
root_P88 <- -1.4264 + 1.2593*root_P50
wb <- hydraulics_psi2Weibull(psi50 = root_P50, psi88 = root_P88)
root_P12 <- medfate::hydraulics_xylemPsi(0.88,1, wb["c"], wb["d"])
x1s$paramsTranspiration$VCroot_P50 <- root_P50
x1s$paramsTranspiration$VCroot_slope <- (88.0 - 12.0)/(abs(root_P88) - abs(root_P12))
x1s$paramsTranspiration$VCroot_c <- wb["c"] 
x1s$paramsTranspiration$VCroot_d <- wb["d"]

medfate:::.updateBelow(x1s)

S1s <- spwb(x1s, pue_meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)
saveRDS(S1s, "Rdata/Puechabon/Real_Puechabon_Sperry_segmented.rds")


# Sperry segmented (TLP) --------------------------------------------------
x1st <- x1s
psi50 = -2.5
psi88 = -1.4264 + 1.2593*psi50
wb <- hydraulics_psi2Weibull(psi50 = psi50, psi88 = psi88)
psi12 <- medfate::hydraulics_xylemPsi(0.88,1, wb["c"], wb["d"])
x1st$paramsTranspiration$VCleaf_P50 <- psi50
x1st$paramsTranspiration$VCleaf_slope <- (88.0 - 12.0)/(abs(psi88) - abs(psi12))
x1st$paramsTranspiration$VCleaf_c <- wb["c"]
x1st$paramsTranspiration$VCleaf_d <- wb["d"]
medfate:::.updateBelow(x1st)

S1st <- spwb(x1st, pue_meteo, 
            latitude = pue_latitude, elevation = pue_elevation, 
            slope = pue_slope, aspect = pue_aspect)
saveRDS(S1st, "Rdata/Puechabon/Real_Puechabon_Sperry_segmented_TLP.rds")


extract_spparams(x1)
extract_spparams(x1s)
extract_spparams(x1st)
extract_spparams(x2b)


# Observed data -----------------------------------------------------------

E_data <- readxl::read_excel("Data/Puechabon/Sapflow_day_SMR_2016-2018.xlsx") |>
  dplyr::select(Date, E_Av.gap)|>
  dplyr::rename(dates = Date, 
                E_T1_168 = E_Av.gap) |>
  dplyr::mutate(dates = as.Date(dates),
                E_T1_168 = E_T1_168/S1$spwbInput$above$LAI_live) |>
  as.data.frame()
row.names(E_data) <- as.character(E_data$dates)

wp_data <- read.table("Data/Puechabon/Water_Potential_MIND_Control-1.csv",sep = ";", 
                      dec=",", header = TRUE) |>
  dplyr::mutate(Date = as.Date(Date, "%d/%m/%Y")) |>
  dplyr::filter(Treatment == "Control") |>
  dplyr::select(-Treatment)  |>
  dplyr::group_by(Date) |>
  dplyr::summarise(PD_T1_168 = mean(Pd_pot, na.rm=TRUE),
            PD_T1_168_err = sqrt(var(Pd_pot, na.rm=TRUE)),
            MD_T1_168 = mean(Md_pot, na.rm=TRUE),
            MD_T1_168_err = sqrt(var(Md_pot, na.rm=TRUE))) |>
  as.data.frame()
row.names(wp_data) <- as.character(wp_data$Date)


library(GA)

# Calibration (non-segmented) ---------------------------------------------
opt_function_ns <- function(par) {
  P50 <- par[1]
  slope <- par[2]
  x1st <- x1 # Non-segmented
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb <- hydraulics_psi2Weibull(psi50 = P50, psi88 = psi88)
  x1st$paramsTranspiration$VCleaf_P50 <- P50
  x1st$paramsTranspiration$VCleaf_slope <- slope
  x1st$paramsTranspiration$VCleaf_c <- wb["c"]
  x1st$paramsTranspiration$VCleaf_d <- wb["d"]
  x1st$paramsTranspiration$VCstem_P50 <- P50
  x1st$paramsTranspiration$VCstem_slope <- slope
  x1st$paramsTranspiration$VCstem_c <- wb["c"]
  x1st$paramsTranspiration$VCstem_d <- wb["d"]
  x1st$paramsTranspiration$VCroot_P50 <- P50
  x1st$paramsTranspiration$VCroot_slope <- slope
  x1st$paramsTranspiration$VCroot_c <- wb["c"]
  x1st$paramsTranspiration$VCroot_d <- wb["d"]
  medfate:::.updateBelow(x1st)
  x1st$control$verbose <- FALSE
  S <- spwb(x1st, pue_meteo, 
            latitude = pue_latitude, elevation = pue_elevation, 
            slope = pue_slope, aspect = pue_aspect)
  mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T1_168", metric = "MAE.rel")
  stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T1_168")
  mae_wp <- mean(stats_wp$MAE.rel)
  mae <- (mae_E + mae_wp)/2
  rm(S)
  cat(paste0("P50 = ", P50, " slope = ", slope, " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  return(-mae)
}
# TEST
opt_function_ns(c(-4.5,21))
g_ns <- ga(type = "real-valued",
        fitness = opt_function_ns,
        lower = c(-5, 10), upper = c(-1,50),
        popSize = 20,
        maxiter = 20,
        optim = FALSE,
        keepBest = TRUE)

# opt = c(-2.488121, 13.36921)
# MAE = -21.35324

P50 <- -2.488121
slope <- 13.36921
x1c <- x1 # Non-segmented
psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb <- hydraulics_psi2Weibull(psi50 = P50, psi88 = psi88)
x1c$paramsTranspiration$VCleaf_P50 <- P50
x1c$paramsTranspiration$VCleaf_slope <- slope
x1c$paramsTranspiration$VCleaf_c <- wb["c"]
x1c$paramsTranspiration$VCleaf_d <- wb["d"]
x1c$paramsTranspiration$VCstem_P50 <- P50
x1c$paramsTranspiration$VCstem_slope <- slope
x1c$paramsTranspiration$VCstem_c <- wb["c"]
x1c$paramsTranspiration$VCstem_d <- wb["d"]
x1c$paramsTranspiration$VCroot_P50 <- P50
x1c$paramsTranspiration$VCroot_slope <- slope
x1c$paramsTranspiration$VCroot_c <- wb["c"]
x1c$paramsTranspiration$VCroot_d <- wb["d"]
medfate:::.updateBelow(x1c)
S1c <- spwb(x1c, pue_meteo, 
          latitude = pue_latitude, elevation = pue_elevation, 
          slope = pue_slope, aspect = pue_aspect)
saveRDS(S1c, "Rdata/Puechabon/Real_Puechabon_Sperry_calibrated.rds")

# Calibration (segmented) -------------------------------------------------------------
opt_function <- function(par) {
  P50 <- par[1]
  slope <- par[2]
  x1st <- x1s
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb <- hydraulics_psi2Weibull(psi50 = P50, psi88 = psi88)
  x1st$paramsTranspiration$VCleaf_P50 <- P50
  x1st$paramsTranspiration$VCleaf_slope <- slope
  x1st$paramsTranspiration$VCleaf_c <- wb["c"]
  x1st$paramsTranspiration$VCleaf_d <- wb["d"]
  medfate:::.updateBelow(x1st)
  x1st$control$verbose <- FALSE
  S <- spwb(x1st, pue_meteo, 
            latitude = pue_latitude, elevation = pue_elevation, 
            slope = pue_slope, aspect = pue_aspect)
  mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T1_168", metric = "MAE.rel")
  stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T1_168")
  mae_wp <- mean(stats_wp$MAE.rel)
  mae <- (mae_E + mae_wp)/2
  rm(S)
  cat(paste0("P50 = ", P50, " slope = ", slope, " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  return(-mae)
}
# TEST
opt_function(c(-4.5,21))

g <- ga(type = "real-valued",
        fitness = opt_function,
        lower = c(-5, 10), upper = c(-1,50),
        popSize = 20,
        maxiter = 20,
        optim = FALSE,
        keepBest = TRUE)

# g <- optim(c(-4.5,20), fn = opt_function, lower = c(-5, 10), upper = c(-1,50), control = list(abstol = 0.0001)) 

# opt = c(-2.514235, 18.19012)
# MAE = -22.24726
x1sc <- x1s
psi50 <- -2.514235 # g$par[1]
slope <- 18.19012 # g$par[2]
psi88 <- psi50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb <- hydraulics_psi2Weibull(psi50 = psi50, psi88 = psi88)
x1sc$paramsTranspiration$VCleaf_P50 <- psi50
x1sc$paramsTranspiration$VCleaf_slope <- slope
x1sc$paramsTranspiration$VCleaf_c <- wb["c"]
x1sc$paramsTranspiration$VCleaf_d <- wb["d"]
medfate:::.updateBelow(x1sc)

S1sc <- spwb(x1sc, pue_meteo, 
             latitude = pue_latitude, elevation = pue_elevation, 
             slope = pue_slope, aspect = pue_aspect)
saveRDS(S1sc, "Rdata/Puechabon/Real_Puechabon_Sperry_segmented_calibrated.rds")

# Plots -------------------------------------------------------------------
library(ggplot2)
library(medfate)
library(cowplot)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry.rds")
S1c <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry_calibrated.rds")
S2b <- readRDS("Rdata/Puechabon/Real_Puechabon_Sureau_Baldocchi.rds")
S1s <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry_segmented.rds")
S1sc <- readRDS("Rdata/Puechabon/Real_Puechabon_Sperry_segmented_calibrated.rds")
p1 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-not-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S1c, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-not-segmented-cal")+theme(legend.position = c(0.8,0.8))
p4 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p5 <- plot(S1sc, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented-cal")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/SoilPsi_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/HydraulicRedistribution_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-8,0))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/LeafPsiRange_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.0))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.0))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/Transpiration_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/StomatalConductance_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/StemPLC_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/LeafPLC_Puechabon_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/SoilPlantConductance_Puechabon_Real.png", p, width = 8, height = 11)

# Evaluation --------------------------------------------------------------
wp_evaluation <- function(S, wp_data, E_data, title) {
  p1 <- evaluation_plot(S, E_data, type="E", cohort = "T1_168")+ ylim(c(0,2))+
    labs(title = title, subtitle = "Sap flux Quercus ilex")+ ylab("Transpiration (kg/m2)")+
    theme_classic()+theme(legend.position = c(0.9,0.85))
  p2<- evaluation_plot(S, wp_data, type="WP", cohort = "T1_168")+ ylim(c(-10,0))+
    labs(subtitle = "Leaf WP Quercus ilex", title = "")+theme_classic()+ theme(legend.position =c(0.1,0.31))
  # p3 <- plot(S, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+
  #   theme_classic()+ theme(legend.position = "none")+labs(title=title, subtitle="PLC")
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1.5,1.5)))  
}
p1 <-wp_evaluation(S2b, wp_data, E_data,"Sureau")
p2 <-wp_evaluation(S1, wp_data, E_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation(S1c, wp_data, E_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation(S1s, wp_data, E_data,"Sperry segmented (measured)")
p5 <-wp_evaluation(S1sc, wp_data, E_data,"Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Puechabon_Real/Evaluation_Puechabon_Real.png", p, width = 17, height = 12, bg = "white")

S1_stats <- evaluation_stats(S1, E_data, type="E", cohort = "T1_168")
S1c_stats <- evaluation_stats(S1c, E_data, type="E", cohort = "T1_168")
S1s_stats <- evaluation_stats(S1s, E_data, type="E", cohort = "T1_168")
S1sc_stats <- evaluation_stats(S1sc, E_data, type="E", cohort = "T1_168")
S2b_stats <- evaluation_stats(S2b, E_data, type="E", cohort = "T1_168")

S2b_stats
S1_stats
S1c_stats
S1s_stats
S1sc_stats
