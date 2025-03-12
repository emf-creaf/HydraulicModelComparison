# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_FontBlanche.R")
source("Rscripts/0_ExtractParams.R")

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
  filter(dates >= as.Date("2016-01-01"),
         dates <= as.Date("2018-12-31"))



# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$stemCavitationRecovery <- "annual"
control$leafCavitationRecovery <- "annual"
control$leafCavitationEffects <- FALSE
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- FALSE
control$sunlitShade <- FALSE

## Total overlap
control$rhizosphereOverlap <- "none"
x1 <- fontblanche_input(control)
S1 <- spwb(x1, fb_meteo,
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1, "Rdata/FontBlanche/Real_FontBlanche_Sperry.rds")


## Partial overlap
# control$rhizosphereOverlap <- "partial"
# x1p <- fontblanche_input(control)
# S1p <- spwb(x1p, fb_meteo, 
#             latitude = fb_latitude, elevation = fb_elevation,
#             slope = fb_slope, aspect = fb_aspect)
# saveRDS(S1p, "Rdata/FontBlanche/Real_FontBlanche_Sperry_partialpools.rds")

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
control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001
control$stomatalSubmodel <- "Baldocchi"
control$rhizosphereOverlap <- "none"
x2b <- fontblanche_input(control)
S2b <- spwb(x2b, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2b, "Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_no_overlap.rds")


## Partial overlap
control$rhizosphereOverlap <- "partial"
x2p <- fontblanche_input(control)
S2p <- spwb(x2p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2p, "Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_partial_overlap.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$rhizosphereOverlap <- "none"
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
# From Limousin et al. (2022)
# Quercus ilex Kleaf P50 = -4.5 P88 = -7.0
# P. halepensis P50 = -2 P88 = -3.2
psi50 <- c(-2.0, -4.5)
slope <- c(40, 20)
psi88 <- psi50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb_1 <- hydraulics_psi2Weibull(psi50 = psi50[1], psi88 = psi88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = psi50[2], psi88 = psi88[2])
psi12 <- c(medfate::hydraulics_xylemPsi(0.88,1, wb_1["c"], wb_1["d"]),
           medfate::hydraulics_xylemPsi(0.88,1, wb_2["c"], wb_2["d"]))
x1s$paramsTranspiration$VCleaf_P50 <- psi50
x1s$paramsTranspiration$VCleaf_slope <- slope
x1s$paramsTranspiration$VCleaf_c <- c(wb_1["c"], wb_2["c"]) 
x1s$paramsTranspiration$VCleaf_d <- c(wb_1["d"], wb_2["d"])
# Bartlett et al. (2016)
# Root P50 = 0.4892 + 0.742 * Stem_P50
root_P50 <- 0.4892 + 0.742*x2b$paramsTranspiration$VCstem_P50
# P88,stem=−1.4264+1.2593* P50,stem
root_P88 <- -1.4264 + 1.2593*root_P50
wb_1 <- hydraulics_psi2Weibull(psi50 = root_P50[1], psi88 = root_P88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = root_P50[2], psi88 = root_P88[2])
root_P12 <- c(medfate::hydraulics_xylemPsi(0.88,1, wb_1["c"], wb_1["d"]),
              medfate::hydraulics_xylemPsi(0.88,1, wb_2["c"], wb_2["d"]))
x1s$paramsTranspiration$VCroot_P50 <- root_P50
x1s$paramsTranspiration$VCroot_slope <- (88.0 - 12.0)/(abs(root_P88) - abs(root_P12))
x1s$paramsTranspiration$VCroot_c <- c(wb_1["c"], wb_2["c"]) 
x1s$paramsTranspiration$VCroot_d <- c(wb_1["d"], wb_2["d"])
medfate:::.updateBelow(x1s)
S1s <- spwb(x1s, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1s, "Rdata/FontBlanche/Real_FontBlanche_Sperry_segmented.rds")

# Sperry segmented (TLP) --------------------------------------------------
x1st <- x1s
psi50 = c(-2.7,-2.5)
psi88 = -1.4264 + 1.2593*psi50
wb_1 <- hydraulics_psi2Weibull(psi50 = psi50[1], psi88 = psi88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = psi50[2], psi88 = psi88[2])
psi12 <- c(medfate::hydraulics_xylemPsi(0.88,1, wb_1["c"], wb_2["d"]),
           medfate::hydraulics_xylemPsi(0.88,1, wb_1["c"], wb_2["d"]))
x1st$paramsTranspiration$VCleaf_P50 <- psi50
x1st$paramsTranspiration$VCleaf_slope <- (88.0 - 12.0)/(abs(psi88) - abs(psi12))
x1st$paramsTranspiration$VCleaf_c <- c(wb_1["c"], wb_2["c"])
x1st$paramsTranspiration$VCleaf_d <- c(wb_1["d"], wb_2["d"])
medfate:::.updateBelow(x1st)

S1st <- spwb(x1st, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1st, "Rdata/FontBlanche/Real_FontBlanche_Sperry_segmented_TLP.rds")



extract_spparams(x1)
extract_spparams(x1s)
extract_spparams(x1st)
extract_spparams(x2b)


# Observed data -----------------------------------------------------------
pot_input <- read.table("Data/FontBlanche/Fontblanche_Potentiels_Hydriques_recapitulatif.txt", header = TRUE)
pot_pine <- pot_input |>
  dplyr::filter(parcelle=="TD", espece == "p") |>
  dplyr::group_by(Date) |>
  dplyr::summarise(PD_T1_148 = mean(Pb), PD_T1_148_err = sqrt(var(Pb)),
                   MD_T1_148 = mean(Pm), MD_T1_148_err = sqrt(var(Pm))) 
pot_oak <- pot_input |>
  dplyr::filter(parcelle=="TD", espece == "c") |>
  dplyr::group_by(Date) |>
  dplyr::summarise(PD_T2_168 = mean(Pb), PD_T2_168_err = sqrt(var(Pb)),
                   MD_T2_168 = mean(Pm), MD_T2_168_err = sqrt(var(Pm))) 
wp_data <- dplyr::left_join(pot_pine, pot_oak, by="Date")|>
  dplyr::mutate(dates = as.Date(Date, format = "%Y-%m-%d"))|>
  as.data.frame()
row.names(wp_data) <- as.character(wp_data$dates)

# pot_20170620 <- read_table("Data/FontBlanche/Fontblanche_Potentiels_Hydriques_20170620.txt", skip = 11) 
# pot_20170712 <- read_table("Data/FontBlanche/Fontblanche_Potentiels_Hydriques_20170712.txt", skip = 11) 
# pot_20170830 <- read_table("Data/FontBlanche/Fontblanche_Potentiels_Hydriques_20170830.txt", skip = 11) 
# 
# extract_wp <- function(x, date) {
#   pot_pine <- x |>
#     dplyr::filter(is.na(remarques), parcelle=="TD") |>
#     dplyr::mutate(sp = substr(arbre,1,1)) |>
#     dplyr::filter(sp == "p") |>
#     dplyr::summarise(PD_T1_148 = mean(Pb), PD_T1_148_err = sqrt(var(Pb)),
#                      MD_T1_148 = mean(Pm), MD_T1_148_err = sqrt(var(Pm))) 
#   pot_oak <- x |>
#     dplyr::filter(is.na(remarques), parcelle=="TD") |>
#     dplyr::mutate(sp = substr(arbre,1,1)) |>
#     dplyr::filter(sp == "c") |>
#     dplyr::summarise(PD_T2_168 = mean(Pb), PD_T2_168_err = sqrt(var(Pb)),
#                      MD_T2_168 = mean(Pm), MD_T2_168_err = sqrt(var(Pm))) 
#   pot <- dplyr::bind_cols(pot_pine, pot_oak)|>
#     dplyr::mutate(dates = as.Date(date))
#   return(pot)
# }
# wp_data <- bind_rows(extract_wp(pot_20170620, "2017-06-20"),
#                      extract_wp(pot_20170712, "2017-07-12"),
#                      extract_wp(pot_20170830, "2017-08-30")) |>
#   as.data.frame()
# row.names(wp_data) <- as.character(wp_data$dates)

FR_FBn_VG_Sap_Flow_3_L2C <- read_delim("Data/FontBlanche/FR-FBn_VG_Sap_Flow_3_L2C.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE,
                                       na = c("-9999"))

E_data <- FR_FBn_VG_Sap_Flow_3_L2C |>
  dplyr::select(TIMESTAMP, SAP_FLOW_TD_P, SAP_FLOW_TD_C) |>
  dplyr::rename(dates = TIMESTAMP,
                E_T1_148 = SAP_FLOW_TD_P, 
                E_T2_168 = SAP_FLOW_TD_C)|>
  dplyr::mutate(dates = as.Date(dates, format = "%d/%m/%Y"))|>
  as.data.frame()
E_data$E_T1_148[E_data$E_T1_148==-9999] <- NA
E_data$E_T2_168[E_data$E_T2_168==-9999] <- NA
lai <- x1$above$LAI_live
E_data$E_T1_148 <- E_data$E_T1_148/lai[1]
E_data$E_T2_168 <- E_data$E_T2_168/lai[2]
row.names(E_data) <- as.character(E_data$dates)


# Calibration (non-segmented) ---------------------------------------------
library(GA)
opt_function_ns <- function(par, cal_species = "Pinus halepensis") {
  if(cal_species =="Pinus halepensis") {
    P50 <- c(par[1], -2.514235)
    slope <- c(par[2], 18.19012)
  } else {
    P50 <- c(-2.507, par[1])
    slope <- c(33.83, par[2])
  }
  x1st <- x1 # Non-segmented
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb_1 <- hydraulics_psi2Weibull(psi50 = P50[1], psi88 = psi88[1])
  wb_2 <- hydraulics_psi2Weibull(psi50 = P50[2], psi88 = psi88[2])
  x1st$paramsTranspiration$VCleaf_P50 <- P50
  x1st$paramsTranspiration$VCleaf_slope <- slope
  x1st$paramsTranspiration$VCleaf_c[1] <- wb_1["c"]
  x1st$paramsTranspiration$VCleaf_d[1] <- wb_1["d"]
  x1st$paramsTranspiration$VCleaf_c[2] <- wb_2["c"]
  x1st$paramsTranspiration$VCleaf_d[2] <- wb_2["d"]
  x1st$paramsTranspiration$VCstem_P50 <- P50
  x1st$paramsTranspiration$VCstem_slope <- slope
  x1st$paramsTranspiration$VCstem_c[1] <- wb_1["c"]
  x1st$paramsTranspiration$VCstem_d[1] <- wb_1["d"]
  x1st$paramsTranspiration$VCstem_c[2] <- wb_2["c"]
  x1st$paramsTranspiration$VCstem_d[2] <- wb_2["d"]
  x1st$paramsTranspiration$VCroot_P50 <- P50
  x1st$paramsTranspiration$VCroot_slope <- slope
  x1st$paramsTranspiration$VCroot_c[1] <- wb_1["c"]
  x1st$paramsTranspiration$VCroot_d[1] <- wb_1["d"]
  x1st$paramsTranspiration$VCroot_c[2] <- wb_2["c"]
  x1st$paramsTranspiration$VCroot_d[2] <- wb_2["d"]
  medfate:::.updateBelow(x1st)
  x1st$control$verbose <- FALSE
  S <- spwb(x1st, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
  if(cal_species == "Pinus halepensis") {
    mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T1_148", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T1_148")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[1], " slope = ", slope[1], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  } else {
    mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T2_168", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T2_168")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[2], " slope = ", slope[2], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  }
  rm(S)
  return(-mae)
}

# Pinus halepensis
# opt_function_ns(c(-4.8,46), cal_species = "Pinus halepensis") # Test
# g_ns_ph <- ga(type = "real-valued",
#               fitness = opt_function_ns,
#               lower = c(-5, 10), upper = c(-1,50),
#               popSize = 20,
#               maxiter = 20,
#               optim = FALSE,
#               keepBest = TRUE,
#               cal_species = "Pinus halepensis")
# opt = c(-2.507, 33.83)
# MAE = 38.70

# Quercus ilex
# opt_function_ns(c(-4.8,46), cal_species = "Quercus ilex") # Test
g_ns_qi <- ga(type = "real-valued",
              fitness = opt_function_ns,
              lower = c(-5, 10), upper = c(-1,50),
              popSize = 20,
              maxiter = 20,
              optim = FALSE,
              keepBest = TRUE,
              cal_species = "Quercus ilex")

# Simulation with calibrated values
P50 <- c(-2.507, -2.514235)
slope <- c(33.83, 18.19012)
x1c <- x1 # Non-segmented
psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb_1 <- hydraulics_psi2Weibull(psi50 = P50[1], psi88 = psi88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = P50[2], psi88 = psi88[2])
x1c$paramsTranspiration$VCleaf_P50 <- P50
x1c$paramsTranspiration$VCleaf_slope <- slope
x1c$paramsTranspiration$VCleaf_c[1] <- wb_1["c"]
x1c$paramsTranspiration$VCleaf_d[1] <- wb_1["d"]
x1c$paramsTranspiration$VCleaf_c[2] <- wb_2["c"]
x1c$paramsTranspiration$VCleaf_d[2] <- wb_2["d"]
x1c$paramsTranspiration$VCstem_P50 <- P50
x1c$paramsTranspiration$VCstem_slope <- slope
x1c$paramsTranspiration$VCstem_c[1] <- wb_1["c"]
x1c$paramsTranspiration$VCstem_d[1] <- wb_1["d"]
x1c$paramsTranspiration$VCstem_c[2] <- wb_2["c"]
x1c$paramsTranspiration$VCstem_d[2] <- wb_2["d"]
x1c$paramsTranspiration$VCroot_P50 <- P50
x1c$paramsTranspiration$VCroot_slope <- slope
x1c$paramsTranspiration$VCroot_c[1] <- wb_1["c"]
x1c$paramsTranspiration$VCroot_d[1] <- wb_1["d"]
x1c$paramsTranspiration$VCroot_c[2] <- wb_2["c"]
x1c$paramsTranspiration$VCroot_d[2] <- wb_2["d"]
medfate:::.updateBelow(x1c)
S1c <- spwb(x1c, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S1c, "Rdata/FontBlanche/Real_FontBlanche_Sperry_calibrated.rds")

# Calibration (segmented) -------------------------------------------------------------
opt_function_s <- function(par, cal_species = "Pinus halepensis") {
  if(cal_species =="Pinus halepensis") {
    P50 <- c(par[1], -2.514235)
    slope <- c(par[2], 18.19012)
  } else {
    P50 <- c(-2.2811, par[1])
    slope <- c(35.690, par[2])
  }
  x1st <- x1s
  psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
  wb_1 <- hydraulics_psi2Weibull(psi50 = P50[1], psi88 = psi88[1])
  wb_2 <- hydraulics_psi2Weibull(psi50 = P50[2], psi88 = psi88[2])
  x1st$paramsTranspiration$VCleaf_P50 <- P50
  x1st$paramsTranspiration$VCleaf_slope <- slope
  x1st$paramsTranspiration$VCleaf_c[1] <- wb_1["c"]
  x1st$paramsTranspiration$VCleaf_d[1] <- wb_1["d"]
  x1st$paramsTranspiration$VCleaf_c[2] <- wb_2["c"]
  x1st$paramsTranspiration$VCleaf_d[2] <- wb_2["d"]
  medfate:::.updateBelow(x1st)
  x1st$control$verbose <- FALSE
  S <- spwb(x1st, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation, 
            slope = fb_slope, aspect = fb_aspect)

  if(cal_species =="Pinus halepensis") {
    mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T1_148", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T1_148")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[1], " slope = ", slope[1], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  } else {
    mae_E <- evaluation_metric(S, E_data, type="E", cohort = "T2_168", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, wp_data, type="WP", cohort = "T2_168")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[2], " slope = ", slope[2], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  }
  rm(S)
  return(-1.0*mae)
}

# Pinus halepensis
opt_function_s(c(-2.0,40), cal_species = "Pinus halepensis") #Test
g_s_ph <- ga(type = "real-valued",
           fitness = opt_function_s,
           lower = c(-5, 10), upper = c(-1,50),
           popSize = 20,
           maxiter = 20,
           optim = FALSE,
           keepBest = TRUE,
           cal_species = "Pinus halepensis")
# opt <- c(-2.28105805384542,  35.689843700499)
# MAE <- 38.65499

# Quercus ilex
opt_function_s(c(-2.0,40), cal_species = "Quercus ilex") #Test
g_s_qi <- ga(type = "real-valued",
             fitness = opt_function_s,
             lower = c(-5, 10), upper = c(-1,50),
             popSize = 20,
             maxiter = 20,
             optim = FALSE,
             keepBest = TRUE,
             cal_species = "Quercus ilex")


# Simulation with calibrated values
P50 <- c(-2.2811, -2.514235)
slope <- c(35.690, 18.19012)
x1sc <- x1s # Segmented
psi88 <- P50  + log((100.0/88.0)-1.0)*(25.0/slope)
wb_1 <- hydraulics_psi2Weibull(psi50 = P50[1], psi88 = psi88[1])
wb_2 <- hydraulics_psi2Weibull(psi50 = P50[2], psi88 = psi88[2])
x1sc$paramsTranspiration$VCleaf_P50 <- P50
x1sc$paramsTranspiration$VCleaf_slope <- slope
x1sc$paramsTranspiration$VCleaf_c[1] <- wb_1["c"]
x1sc$paramsTranspiration$VCleaf_d[1] <- wb_1["d"]
x1sc$paramsTranspiration$VCleaf_c[2] <- wb_2["c"]
x1sc$paramsTranspiration$VCleaf_d[2] <- wb_2["d"]
medfate:::.updateBelow(x1sc)

S1sc <- spwb(x1sc, fb_meteo, 
             latitude = fb_latitude, elevation = fb_elevation, 
             slope = fb_slope, aspect = fb_aspect)
saveRDS(S1sc, "Rdata/FontBlanche/Real_FontBlanche_Sperry_segmented_calibrated.rds")

# Plots -------------------------------------------------------------------
library(medfate)
library(cowplot)
library(ggplot2)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/FontBlanche/Real_FontBlanche_Sperry.rds")
S1c <- readRDS("Rdata/FontBlanche/Real_FontBlanche_Sperry_calibrated.rds")
S1s <- readRDS("Rdata/FontBlanche/Real_FontBlanche_Sperry_segmented.rds")
S1sc <- readRDS("Rdata/FontBlanche/Real_FontBlanche_Sperry_segmented_calibrated.rds")
S2b <- readRDS("Rdata/FontBlanche/Real_FontBlanche_Sureau_Baldocchi_no_overlap.rds")

p1 <- plot(S2b, "SoilPsi")+ylim(c(-7,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-not-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S1c, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-not-segmented-cal")+theme(legend.position = c(0.8,0.8))
p4 <- plot(S1s, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p5 <- plot(S1sc, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-segmented-cal")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/SoilPsi_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/HydraulicRedistribution_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sureau")
p2 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/LeafPsiRange_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sureau")
p2 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/Transpiration_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sureau")
p2 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/StomatalConductance_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sureau")
p2 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/StemPLC_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p2 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/LeafPLC_FontBlanche_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/SoilPlantConductance_FontBlanche_Real.png", p, width = 8, height = 11)


wp_evaluation_oak <- function(S, wp_data, E_data, title) {
  p1 <- evaluation_plot(S, E_data, type="E", cohort = "T2_168")+ ylim(c(0,2))+
    labs(title = title, subtitle = "Sap flux Quercus ilex")+ ylab("")+
    theme_classic()+theme(legend.position = c(0.9,0.8))
  p2<- evaluation_plot(S, wp_data, type="WP", cohort = "T2_168")+ ylim(c(-10,0))+ylab("")+
    labs(subtitle = "Leaf WP Quercus ilex", title = title)+theme_classic()+ theme(legend.position = c(0.1,0.35))
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1,1)))  
}
p1 <-wp_evaluation_oak(S2b, wp_data, E_data, "Sureau")
p2 <-wp_evaluation_oak(S1, wp_data, E_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation_oak(S1c, wp_data, E_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation_oak(S1s, wp_data, E_data, "Sperry segmented (measured)")
p5 <-wp_evaluation_oak(S1sc, wp_data, E_data, "Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/Evaluation_FontBlanche_Real_Oak.png", p, width = 17, height = 12, bg = "white")

evaluation_stats(S2b, E_data, type="E", cohort = "T2_168")
evaluation_stats(S2b, wp_data, type="WP", cohort = "T2_168")
evaluation_stats(S1, E_data, type="E", cohort = "T2_168")
evaluation_stats(S1, wp_data, type="WP", cohort = "T2_168")
evaluation_stats(S1c, E_data, type="E", cohort = "T2_168")
evaluation_stats(S1c, wp_data, type="WP", cohort = "T2_168")
evaluation_stats(S1s, E_data, type="E", cohort = "T2_168")
evaluation_stats(S1s, wp_data, type="WP", cohort = "T2_168")
evaluation_stats(S1sc, E_data, type="E", cohort = "T2_168")
evaluation_stats(S1sc, wp_data, type="WP", cohort = "T2_168")

wp_evaluation_pine <- function(S, wp_data, E_data, title) {
  p1 <- evaluation_plot(S, E_data, type="E", cohort = "T1_148")+ ylim(c(0,1.5))+ ylab("Transpiration (kg/m2)")+
    labs(title = title, subtitle = "Sap flux Pinus halepensis")+
    theme_classic()+theme(legend.position = c(0.9,0.8))
  p2<- evaluation_plot(S, wp_data, type="WP", cohort = "T1_148")+ylim(c(-8,0))+
    theme_classic()+theme(legend.position = c(0.1,0.35)) + labs(title = title, subtitle = "Leaf WP Pinus halepensis")
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1,1)))  
}
p1 <-wp_evaluation_pine(S2b, wp_data, E_data, "Sureau")
p2 <-wp_evaluation_pine(S1, wp_data, E_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation_pine(S1c, wp_data, E_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation_pine(S1s, wp_data, E_data, "Sperry segmented (measured)")
p5 <-wp_evaluation_pine(S1sc, wp_data, E_data, "Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/FontBlanche_Real/Evaluation_FontBlanche_Real_Pine.png", p, width = 17, height = 12, bg = "white")

evaluation_stats(S2b, E_data, type="E", cohort = "T1_148")
evaluation_stats(S2b, wp_data, type="WP", cohort = "T1_148")
evaluation_stats(S1, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1, wp_data, type="WP", cohort = "T1_148")
evaluation_stats(S1c, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1c, wp_data, type="WP", cohort = "T1_148")
evaluation_stats(S1s, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1s, wp_data, type="WP", cohort = "T1_148")
evaluation_stats(S1sc, E_data, type="E", cohort = "T1_148")
evaluation_stats(S1sc, wp_data, type="WP", cohort = "T1_148")


