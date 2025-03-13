# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0c_Ancillary_Yatir.R")
source("Rscripts/0_ExtractParams.R")

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
  filter(dates >= as.Date("2012-01-01"),
         dates <= as.Date("2014-12-31"))



# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$stemCavitationRecovery <- "annual"
control$leafCavitationRecovery <- "annual"
control$leafCavitationEffects <- TRUE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE
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
leaf_P50 <- -2.0
leaf_slope <- 40
leaf_P88 <- leaf_P50 + log((1/0.88)-1.0)*(25.0/leaf_slope)
wb <- hydraulics_psi2Weibull(psi50 = leaf_P50, psi88 = leaf_P88)
x1s$paramsTranspiration$VCleaf_P50 <- leaf_P50
x1s$paramsTranspiration$VCleaf_slope <- leaf_slope
x1s$paramsTranspiration$VCleaf_c <- wb["c"]
x1s$paramsTranspiration$VCleaf_d <- wb["d"]
# Bartlett et al. (2016)
# Root P50 = 0.4892 + 0.742 * Stem_P50
root_P50 <- 0.4892 + 0.742*x2b$paramsTranspiration$VCstem_P50
# P88,stem=−1.4264+1.2593* P50,stem
root_P88 <- -1.4264 + 1.2593*root_P50
wb <- hydraulics_psi2Weibull(psi50 = root_P50, psi88 = root_P88)
root_P12 <- c(medfate::hydraulics_xylemPsi(0.88,1, wb["c"], wb["d"]))
x1s$paramsTranspiration$VCroot_P50 <- root_P50
x1s$paramsTranspiration$VCroot_slope <- (88.0 - 12.0)/(abs(root_P88) - abs(root_P12))
x1s$paramsTranspiration$VCroot_d <- wb["d"]
x1s$paramsTranspiration$VCroot_c <- wb["c"] 
x1s$paramsTranspiration$VCroot_d <- wb["d"]

S1s <- spwb(x1s, yat_meteo, 
           latitude = yat_latitude, elevation = yat_elevation, 
           slope = yat_slope, aspect = yat_aspect)
saveRDS(S1s, "Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")


# Observed data -----------------------------------------------------------
# sapflow data, está en dm3/h, y el timestep es 30 minutos, así que si dividimos
# entre dos ya tenemos los L en esa media hora. sumamos todo el día y 
# promediamos entre arboles. a continuación multiplicamos por la densidad de la parcela
# y dividimos por el LAI para obtener el valor por superficie de hoja
sapf_data <- read.csv('Data/Yatir/ISR_YAT_YAT_sapf_data.csv')
stand_md <- read.csv('Data/Yatir/ISR_YAT_YAT_stand_md.csv')
transp_data_temp <- sapf_data |>
  dplyr::mutate_at(dplyr::vars(dplyr::starts_with('ISR_YAT_YAT')),
                   dplyr::funs(./2)) |>
  dplyr::mutate(dates = date(as_datetime(TIMESTAMP, tz = 'Europe/Madrid'))) |>
  # sum days
  dplyr::group_by(dates) |>
  dplyr::summarise_at(dplyr::vars(dplyr::starts_with('ISR_YAT_YAT')),
                      dplyr::funs(sum(., na.rm = TRUE))) |>
  # remove the zeroes generated previously
  dplyr::mutate_at(dplyr::vars(dplyr::starts_with('ISR_YAT_YAT')),
                   dplyr::funs(replace(., . == 0, NA)))
transp_data_temp$E_Ph <- rowMeans(transp_data_temp[,2:25], na.rm=TRUE)
transp_data_temp2 <- transp_data_temp |>
  dplyr::select(dates, E_Ph) |>
  dplyr::mutate(E_Ph = E_Ph*(stand_md$st_density/(2.0*10000)))
names(transp_data_temp2)[2] <- paste0("E_", "T1_148")
E_data <- transp_data_temp2 |>
  dplyr::filter(!is.na(dates) & !is.na(E_T1_148)) |>
  as.data.frame()
row.names(E_data) <- as.character(E_data$dates)

ggplot(E_data) +
  geom_line(aes(x=dates, y=E_T1_148))

# Calibration (non-segmented) ---------------------------------------------
library(GA)
opt_function_ns <- function(par) {
  P50 <- c(par[1])
  slope <- c(par[2])
  x1st <- x1 # Non-segmented
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb <- hydraulics_psi2Weibull(psi50 = P50[1], psi88 = psi88[1])
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
  S <- spwb(x1st, yat_meteo,
            latitude = yat_latitude, elevation = yat_elevation,
            slope = yat_slope, aspect = yat_aspect)
  mae <- evaluation_metric(S, E_data, type="E", cohort = "T1_148", metric = "MAE.rel")
  cat(paste0("P50 = ", P50, " slope = ", slope, " MAE = ", mae, "\n"))
  return(-mae)
}

# TEST
opt_function_ns(c(-4.8,46))
g_ns <- ga(type = "real-valued",
           fitness = opt_function_ns,
           lower = c(-5, 10), upper = c(-1,50),
           popSize = 40,
           maxiter = 30,
           optim = FALSE,
           keepBest = TRUE)
saveRDS(g_ns, "Rdata/Yatir/g_ns_ph.rds")

g_ns <- readRDS("Rdata/Yatir/g_ns_ph.rds")
# opt = c(-2.340047, 20.47438)
# MAE = 47.99252

P50 <- g_ns@solution[1]
slope <- g_ns@solution[2]
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
S1c <- spwb(x1c, yat_meteo, 
            latitude = yat_latitude, elevation = yat_elevation, 
            slope = yat_slope, aspect = yat_aspect)
saveRDS(S1c, "Rdata/Yatir/Real_Yatir_Sperry_calibrated.rds")

# Calibration (segmented) -------------------------------------------------------------
library(GA)
opt_function_s <- function(par) {
  P50 <- c(par[1])
  slope <- c(par[2])
  x1st <- x1s
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb <- hydraulics_psi2Weibull(psi50 = P50, psi88 = psi88)
  x1st$paramsTranspiration$VCleaf_P50 <- P50
  x1st$paramsTranspiration$VCleaf_slope <- slope
  x1st$paramsTranspiration$VCleaf_c <- wb["c"]
  x1st$paramsTranspiration$VCleaf_d <- wb["d"]
  medfate:::.updateBelow(x1st)
  x1st$control$verbose <- FALSE
  S <- spwb(x1st, yat_meteo,
            latitude = yat_latitude, elevation = yat_elevation,
            slope = yat_slope, aspect = yat_aspect)
  mae <- evaluation_metric(S, E_data, type="E", cohort = "T1_148", metric = "MAE.rel")
  cat(paste0("P50 = ", P50, " slope = ", slope, " MAE = ", mae, "\n"))
  return(-1.0*mae)
}

opt_function_s(c(-2.0,40))

g_s <- ga(type = "real-valued",
        fitness = opt_function_s,
        lower = c(-5, 10), upper = c(-1,50),
        popSize = 40,
        maxiter = 30,
        optim = FALSE,
        keepBest = TRUE)
saveRDS(g_s, "Rdata/Yatir/g_s_ph.rds")

g_s <- readRDS("Rdata/Yatir/g_s_ph.rds")
# opt <- c(-2.676182,  29.65367)
# MAE <- 48.34733
P50 <- g_s@solution[1]
slope <- g_s@solution[2]
x1sc <- x1s
psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb <- hydraulics_psi2Weibull(psi50 = P50, psi88 = psi88)
x1sc$paramsTranspiration$VCleaf_P50 <- P50
x1sc$paramsTranspiration$VCleaf_slope <- slope
x1sc$paramsTranspiration$VCleaf_c <- wb["c"]
x1sc$paramsTranspiration$VCleaf_d <- wb["d"]
medfate:::.updateBelow(x1sc)

S1sc <- spwb(x1sc, yat_meteo, 
            latitude = yat_latitude, elevation = yat_elevation, 
            slope = yat_slope, aspect = yat_aspect)
saveRDS(S1sc, "Rdata/Yatir/Real_Yatir_Sperry_segmented_calibrated.rds")

# Plots -------------------------------------------------------------------
library(ggplot2)
library(cowplot)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/Yatir/Real_Yatir_Sperry.rds")
S1c <- readRDS("Rdata/Yatir/Real_Yatir_Sperry_calibrated.rds")
S1s <- readRDS("Rdata/Yatir/Real_Yatir_Sperry_segmented.rds")
S1sc <- readRDS("Rdata/Yatir/Real_Yatir_Sperry_segmented_calibrated.rds")
S2b <- readRDS("Rdata/Yatir/Real_Yatir_Sureau_Baldocchi.rds")

p1 <- plot(S2b, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-not-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S1c, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-not-segmented-calivrated")+theme(legend.position = c(0.8,0.8))
p4 <- plot(S1s, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p5 <- plot(S1sc, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry-segmented-cal")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/SoilPsi_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/HydraulicRedistribution_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p4 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-7,0))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4,  p5, nrow = 5)
ggsave2("Plots/Yatir_Real/LeafPsiRange_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/Transpiration_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/StomatalConductance_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/StemPLC_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/LeafPLC_Yatir_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/SoilPlantConductance_Yatir_Real.png", p, width = 8, height = 11)


wp_evaluation <- function(S, E_data, title) {
  p1 <- evaluation_plot(S, E_data, type="E", cohort = "T1_148")+ ylim(c(0,2))+
    labs(title = title, subtitle = "Sap flux Pinus halepensis")+ ylab("")+
    theme_classic()+theme(legend.position = c(0.9,0.8))
  p2<- plot(S, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-6.0,0))+
    labs(subtitle = "Leaf WP Pinus halepensis", title = "")+
    theme_classic()+ theme(legend.position = "none")
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1,1)))  
}
# wp_evaluation <- function(S,title) {
#   p1 <- plot(S, "Transpiration", bySpecies = TRUE)+ylim(c(0,1.0))+
#     labs(title = title, subtitle = "Transpiration Pinus halepensis")+
#     theme_classic()+theme(legend.position = "none")
#   p3 <- plot(S, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+
#     theme_classic()+ theme(legend.position = "none")+
#     labs(title=title, subtitle="PLC")
#   return(cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(2,2,1.3)))  
# }
p1 <-wp_evaluation(S2b, E_data, "Sureau")
p2 <-wp_evaluation(S1, E_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation(S1c, E_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation(S1s, E_data, "Sperry segmented (measured)")
p5 <-wp_evaluation(S1sc, E_data, "Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Yatir_Real/Evaluation_Yatir_Real.png", p, width = 17, height = 12, bg = "white")


evaluation_stats(S2b, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1c, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1s, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1sc, E_data, type="E", cohort = "T1_148")


extract_spparams(x1)
extract_spparams(x1s)
extract_spparams(x2b)

