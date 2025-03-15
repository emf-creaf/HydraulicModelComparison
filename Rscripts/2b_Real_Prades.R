# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_Prades.R")
source("Rscripts/0_ExtractParams.R")

# Terrain -----------------------------------------------------------------
pr_latitude <- 41.33
pr_elevation <- 1018
pr_slope <- 35
pr_aspect <- 8.5

# Soil --------------------------------------------------------------------
pr_soil <- prades_soil()
sum(soil_waterFC(pr_soil, "VG"))
sum(soil_waterExtractable(pr_soil, "VG", -1.5))

# Meteo -------------------------------------------------------------------
pr_meteo <- read.table("Data/Prades/PRADES_meteoData.txt", sep="\t", header = TRUE) |>
  dplyr::mutate(dates = as.Date(dates))|>
  as.data.frame()

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
# control$sunlitShade <- FALSE
control$soilDomains <- "dual"
control$segmentedXylemVulnerability <- FALSE

## Total overlap
control$rhizosphereOverlap <- "none"
x1 <- prades_input(control)
S1 <- spwb(x1, pr_meteo,
            latitude = pr_latitude, elevation = pr_elevation,
            slope = pr_slope, aspect = pr_aspect)
saveRDS(S1, "Rdata/Prades/Real_Prades_Sperry.rds")


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
# control$sunlitShade <- FALSE
control$gs_NightFrac <- 0.001
control$stomatalSubmodel <- "Baldocchi"
control$rhizosphereOverlap <- "none"
x2b <- prades_input(control)
S2b <- spwb(x2b, pr_meteo, 
            latitude = pr_latitude, elevation = pr_elevation,
            slope = pr_slope, aspect = pr_aspect)
saveRDS(S2b, "Rdata/Prades/Real_Prades_Sureau_Baldocchi_no_overlap.rds")

# Sperry (segmented) ------------------------------------------------------
x1s <- x1
x1s$control$rhizosphereOverlap <- "none"
x1s$control$leafCavitationRecovery <- "total"
x1s$control$leafCavitationEffects <- FALSE
# Quercus ilex Kleaf P50 = -4.5 P88 = -7.0
# P. sylvestris P50 = -1.30 P88 = -3.2
psi50 <- c(-1.30, -4.5)
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
S1s <- spwb(x1s, pr_meteo, 
            latitude = pr_latitude, elevation = pr_elevation, 
            slope = pr_slope, aspect = pr_aspect)
saveRDS(S1s, "Rdata/Prades/Real_Prades_Sperry_segmented.rds")



# Observed data -----------------------------------------------------------
measured_data <- read.table("Data/Prades/PRADES_measuredData.txt", sep="\t", header = TRUE)

# Calibration (non-segmented) ---------------------------------------------
library(GA)
opt_function_ns <- function(par, cal_species = "Pinus sylvestris") {
  if(cal_species =="Pinus sylvestris") {
    P50 <- c(par[1], -2.514235)
    slope <- c(par[2], 18.19012)
  } else {
    P50 <- c(-2.73, par[1])
    slope <- c(43.24, par[2])
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
  S <- spwb(x1st, pr_meteo,
            latitude = pr_latitude, elevation = pr_elevation,
            slope = pr_slope, aspect = pr_aspect)
  if(cal_species == "Pinus sylvestris") {
    mae_E <- evaluation_metric(S, measured_data, type="E", cohort = "T1_361", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, measured_data, type="WP", cohort = "T1_361")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[1], " slope = ", slope[1], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  } else {
    mae_E <- evaluation_metric(S, measured_data, type="E", cohort = "T2_394", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, measured_data, type="WP", cohort = "T2_394")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[2], " slope = ", slope[2], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  }
  rm(S)
  return(-mae)
}

# Pinus sylvestris
opt_function_ns(c(-4.8,46), cal_species = "Pinus sylvestris") # Test
g_ns_ps <- ga(type = "real-valued",
              fitness = opt_function_ns,
              lower = c(-5, 10), upper = c(-1,50),
              popSize = 40,
              maxiter = 20,
              optim = FALSE,
              keepBest = TRUE,
              cal_species = "Pinus sylvestris")
saveRDS(g_ns_ps, "Rdata/Prades/g_ns_ps.rds")

g_ns_ps <- readRDS("Rdata/Prades/g_ns_ps.rds")
# opt = c(-2.741812, 40.55099)
# MAE = 40.25751

# Quercus ilex
opt_function_ns(c(-4.8,46), cal_species = "Quercus ilex") # Test
g_ns_qi <- ga(type = "real-valued",
              fitness = opt_function_ns,
              lower = c(-5, 10), upper = c(-1,50),
              popSize = 40,
              maxiter = 20,
              optim = FALSE,
              keepBest = TRUE,
              cal_species = "Quercus ilex")
saveRDS(g_ns_qi, "Rdata/Prades/g_ns_qi.rds")

g_ns_qi <- readRDS("Rdata/Prades/g_ns_qi.rds")
# opt = c(-2.168199, 18.06897)
# MAE = 30.69368

# Simulation with calibrated values
P50 <- c(g_ns_ps@solution[1], g_ns_qi@solution[1])
slope <- c(g_ns_ps@solution[2], g_ns_qi@solution[2])
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
S1c <- spwb(x1c, pr_meteo, 
            latitude = pr_latitude, elevation = pr_elevation, 
            slope = pr_slope, aspect = pr_aspect)
saveRDS(S1c, "Rdata/Prades/Real_Prades_Sperry_calibrated.rds")

# Calibration (segmented) -------------------------------------------------------------
library(GA)
opt_function_s <- function(par, cal_species = "Pinus sylvestris") {
  if(cal_species =="Pinus sylvestris") {
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
  S <- spwb(x1st, pr_meteo,
            latitude = pr_latitude, elevation = pr_elevation,
            slope = pr_slope, aspect = pr_aspect)

  if(cal_species =="Pinus sylvestris") {
    mae_E <- evaluation_metric(S, measured_data, type="E", cohort = "T1_361", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, measured_data, type="WP", cohort = "T1_361")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[1], " slope = ", slope[1], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  } else {
    mae_E <- evaluation_metric(S, measured_data, type="E", cohort = "T2_394", metric = "MAE.rel")
    stats_wp <- evaluation_stats(S, measured_data, type="WP", cohort = "T2_394")
    mae_wp <- mean(stats_wp$MAE.rel)
    mae <- (mae_E + mae_wp)/2
    cat(paste0("P50 = ", P50[2], " slope = ", slope[2], " MAE(E) = ", mae_E, " MAE(wp) = ", mae_wp, " MAE = ", mae, "\n"))
  }
  rm(S)
  return(-1.0*mae)
}

# Pinus sylvestris
opt_function_s(c(-2.0,40), cal_species = "Pinus sylvestris") #Test
g_s_ps <- ga(type = "real-valued",
           fitness = opt_function_s,
           lower = c(-5, 10), upper = c(-1,50),
           popSize = 40,
           maxiter = 20,
           optim = FALSE,
           keepBest = TRUE,
           cal_species = "Pinus sylvestris")
saveRDS(g_s_ps, "Rdata/Prades/g_s_ps.rds")

g_s_ps <- readRDS("Rdata/Prades/g_s_ps.rds")
# opt <- c(-2.559114,  41.60563)
# MAE <- 40.42071

# Quercus ilex
opt_function_s(c(-2.0,40), cal_species = "Quercus ilex") #Test
g_s_qi <- ga(type = "real-valued",
             fitness = opt_function_s,
             lower = c(-5, 10), upper = c(-1,50),
             popSize = 40,
             maxiter = 20,
             optim = FALSE,
             keepBest = TRUE,
             cal_species = "Quercus ilex")
saveRDS(g_s_qi, "Rdata/Prades/g_s_qi.rds")

g_s_qi <- readRDS("Rdata/Prades/g_s_qi.rds")
# opt <- c(-2.015863 ,  23.38618)
# MAE <- 31.52489

# Simulation with calibrated values
P50 <- c(g_s_ps@solution[1], g_s_qi@solution[1])
slope <- c(g_s_ps@solution[2], g_s_qi@solution[2])
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

S1sc <- spwb(x1sc, pr_meteo, 
             latitude = pr_latitude, elevation = pr_elevation, 
             slope = pr_slope, aspect = pr_aspect)
saveRDS(S1sc, "Rdata/Prades/Real_Prades_Sperry_segmented_calibrated.rds")

# Plots -------------------------------------------------------------------
library(medfate)
library(cowplot)
library(ggplot2)
Sys.setlocale("LC_ALL", "en_US.utf8")
S1 <- readRDS("Rdata/Prades/Real_Prades_Sperry.rds")
S1c <- readRDS("Rdata/Prades/Real_Prades_Sperry_calibrated.rds")
S1s <- readRDS("Rdata/Prades/Real_Prades_Sperry_segmented.rds")
S1sc <- readRDS("Rdata/Prades/Real_Prades_Sperry_segmented_calibrated.rds")
S2b <- readRDS("Rdata/Prades/Real_Prades_Sureau_Baldocchi_no_overlap.rds")

p1 <- plot(S2b, "SoilPsi")+ylim(c(-7,0))+labs(title="Sureau")+theme(legend.position = c(0.8,0.8))
p2 <- plot(S1, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-not-segmented")+theme(legend.position = c(0.8,0.8))
p3 <- plot(S1c, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-not-segmented-cal")+theme(legend.position = c(0.8,0.8))
p4 <- plot(S1s, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-segmented")+theme(legend.position = c(0.8,0.8))
p5 <- plot(S1sc, "SoilPsi")+ylim(c(-7,0))+labs(title="Sperry-segmented-cal")+theme(legend.position = c(0.8,0.8))
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/SoilPsi_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "HydraulicRedistribution")+ylim(c(0,0.5))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/HydraulicRedistribution_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sureau")
p2 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-10,0))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/LeafPsiRange_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sureau")
p2 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = c(0.1,0.2))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/Transpiration_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sureau")
p2 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "GSWMax_SL", bySpecies = TRUE)+theme(legend.position = c(0.1,0.2))+ylim(c(0,0.3))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/StomatalConductance_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sureau")+theme_classic()
p2 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-not-segmented")+theme_classic()
p3 <- plot(S1c, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-not-segmented-cal")+theme_classic()
p4 <- plot(S1s, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-segmented")+theme_classic()
p5 <- plot(S1sc, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = c(0.15,0.8))+labs(title="Sperry-segmented-cal")+theme_classic()
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/StemPLC_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+labs(title="Sureau", subtitle = "Stem PLC")+ylab("Stem PLC (%)") + theme_classic()+theme(legend.position = c(0.1,0.8))
p2 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+labs(title="Sperry non-segmented (measured)", subtitle = "Stem PLC")+ylab("Stem PLC (%)")+theme_classic()+theme(legend.position = c(0.1,0.8))
p3 <- plot(S1c, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+labs(title="Sperry non-segmented (calibrated)", subtitle = "Stem PLC")+ylab("Stem PLC (%)")+theme_classic()+theme(legend.position = c(0.1,0.8))
p4 <- plot(S1s, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+labs(title="Sperry segmented (measured)", subtitle = "Stem PLC")+ylab("Stem PLC (%)")+theme_classic()+theme(legend.position = c(0.1,0.8))
p5 <- plot(S1sc, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+labs(title="Sperry segmented (calibrated)", subtitle = "Stem PLC")+ylab("Stem PLC (%)")+theme_classic()+theme(legend.position = c(0.1,0.8))
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/StemPLC_Prades_Real.png", p, width = 8, height = 11)

p1 <- plot(S2b, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sureau")
p2 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented")
p3 <- plot(S1c, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-not-segmented-cal")
p4 <- plot(S1s, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented")
p5 <- plot(S1sc, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = c(0.8,0.8))+labs(title="Sperry-segmented-cal")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/SoilPlantConductance_Prades_Real.png", p, width = 8, height = 11)


wp_evaluation_oak <- function(S, measured_data, title) {
  p1 <- evaluation_plot(S, measured_data, type="E", cohort = "T2_394")+ ylim(c(0,2))+
    labs(title = title, subtitle = "Sap flux Quercus ilex")+ ylab("")+
    theme_classic()+theme(legend.position = c(0.9,0.8))
  p2<- evaluation_plot(S, measured_data, type="WP", cohort = "T2_394")+ ylim(c(-10,0))+ylab("")+
    labs(subtitle = "Leaf WP Quercus ilex", title = title)+theme_classic()+ theme(legend.position = c(0.1,0.35))
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1,1)))
}
p1 <-wp_evaluation_oak(S2b, measured_data, "Sureau")
p2 <-wp_evaluation_oak(S1, measured_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation_oak(S1c, measured_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation_oak(S1s, measured_data, "Sperry segmented (measured)")
p5 <-wp_evaluation_oak(S1sc, measured_data, "Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/Evaluation_Prades_Real_Oak.png", p, width = 17, height = 12, bg = "white")

evaluation_stats(S2b, measured_data, type="E", cohort = "T2_394")
evaluation_stats(S2b, measured_data, type="WP", cohort = "T2_394")
evaluation_stats(S1, measured_data, type="E", cohort = "T2_394")
evaluation_stats(S1, measured_data, type="WP", cohort = "T2_394")
evaluation_stats(S1c, measured_data, type="E", cohort = "T2_394")
evaluation_stats(S1c, measured_data, type="WP", cohort = "T2_394")
evaluation_stats(S1s, measured_data, type="E", cohort = "T2_394")
evaluation_stats(S1s, measured_data, type="WP", cohort = "T2_394")
evaluation_stats(S1sc, measured_data, type="E", cohort = "T2_394")
evaluation_stats(S1sc, measured_data, type="WP", cohort = "T2_394")

wp_evaluation_pine <- function(S, measured_data, title) {
  p1 <- evaluation_plot(S, measured_data, type="E", cohort = "T1_361")+ ylim(c(0,1.5))+ ylab("Transpiration (kg/m2)")+
    labs(title = title, subtitle = "Sap flux Pinus sylvestris")+
    theme_classic()+theme(legend.position = c(0.9,0.8))
  p2<- evaluation_plot(S, measured_data, type="WP", cohort = "T1_361")+ylim(c(-6,0))+
    theme_classic()+theme(legend.position = c(0.07,0.35)) + labs(title = title, subtitle = "Leaf WP Pinus sylvestris")
  return(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1,1)))
}
p1 <-wp_evaluation_pine(S2b, measured_data, "Sureau")
p2 <-wp_evaluation_pine(S1, measured_data, "Sperry non-segmented (measured)")
p3 <-wp_evaluation_pine(S1c, measured_data, "Sperry non-segmented (calibrated)")
p4 <-wp_evaluation_pine(S1s, measured_data, "Sperry segmented (measured)")
p5 <-wp_evaluation_pine(S1sc, measured_data, "Sperry segmented (calibrated)")
p <-plot_grid(p1, p2, p3, p4, p5, nrow = 5)
ggsave2("Plots/Prades_Real/Evaluation_Prades_Real_Pine.png", p, width = 17, height = 12, bg = "white")

evaluation_stats(S2b, measured_data, type="E", cohort = "T1_361")
evaluation_stats(S2b, measured_data, type="WP", cohort = "T1_361")
evaluation_stats(S1, measured_data, type="E", cohort = "T1_361")
evaluation_stats(S1, measured_data, type="WP", cohort = "T1_361")
evaluation_stats(S1c, measured_data, type="E", cohort = "T1_361")
evaluation_stats(S1c, measured_data, type="WP", cohort = "T1_361")
evaluation_stats(S1s, measured_data, type="E", cohort = "T1_361")
evaluation_stats(S1s, measured_data, type="WP", cohort = "T1_361")
evaluation_stats(S1sc, measured_data, type="E", cohort = "T1_361")
evaluation_stats(S1sc, measured_data, type="WP", cohort = "T1_361")


extract_spparams(x1)
extract_spparams(x1s)
extract_spparams(x1sc)
extract_spparams(x2b)
