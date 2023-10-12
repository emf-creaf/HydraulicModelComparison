# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)

FrFBn_Qilex_LAI2 <- read_table("Data/Puechabon/FrFBn_Qilex_LAI2.csv")

k_plant = 0.8
k_RSApo1 <- 0.364
k_RSApo2 <- 0.903
k_RSApo3 <- 0.732
# k_RSApo <- k_RSApo1 + k_RSApo2 + k_RSApo3
k_RSApo <- 2
k_SLApo <- 4
k_LSym <- 2
1/(1/k_RSApo + 1/k_SLApo + 1/k_LSym)

# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270

# Soil --------------------------------------------------------------------
pue_soil_df <- defaultSoilParams(n = 3)
pue_soil_df$widths <- c(200, 800, 3000)
pue_soil_df$rfc <- c(75, 82, 94)
pue_soil_df$clay <- rep(39, 3)
pue_soil_df$sand <- rep(26, 3)
pue_soil_df$om <- c(6,3,1)
pue_soil_df$bd <-rep(1.45, 3)

pue_soil <- soil(pue_soil_df)
pue_soil$VG_theta_sat <- rep(0.5, 3)
pue_soil$VG_theta_res <- rep(0.098, 3)
pue_soil$VG_n <- rep(1.55, 3)
# 1 cm = 0.00009804139432 MPa
pue_soil$VG_alpha <- rep(0.0005/0.00009804139432, 3)
pue_soil$Ksat <- rep(655.2934*1.69, 3)

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


# Forest ------------------------------------------------------------------

pue_forest <- emptyforest(1)
pue_forest$treeData$Species = "Quercus ilex"
pue_forest$treeData$DBH = 9.11
pue_forest$treeData$Height = 530.2
pue_forest$treeData$N = 1750
pue_forest$treeData$Z50 = 306.5
pue_forest$treeData$Z95 = 1000
pue_forest$treeData$LAI = 2.2


# Sperry initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- FALSE
control$cavitationRefill <- "annual"
control$bareSoilEvaporation <- TRUE
control$sapFluidityVariation <- FALSE


#Initialize input
x1 <- forest2spwbInput(pue_forest, pue_soil, SpParamsMED, control)

x1$paramsInterception$kDIR <- 0.4
# K	0.5
x1$paramsInterception$kPAR <- 0.5
# canopyStorageParam	1.5
x1$paramsInterception$g <- 1.5/2.2
# P50_VC_Leaf	-6.4
# slope_VC_Leaf	30
P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
wb <- hydraulics_psi2Weibull(-6.4, P88)
# wbl <- hydraulics_psi2Weibull(-4.5, -5.5)
wbl <- hydraulics_psi2Weibull(-3.5, -4.5)
x1$paramsTranspiration$VCleaf_c <- wbl["c"]
x1$paramsTranspiration$VCstem_c <- wb["c"]
x1$paramsTranspiration$VCroot_c <- wb["c"]
x1$paramsTranspiration$VCleaf_d <- wbl["d"]
x1$paramsTranspiration$VCstem_d <- wb["d"]
x1$paramsTranspiration$VCroot_d <- wb["d"]
# gmin20	2
# gmin_S	2
x1$paramsTranspiration$Gswmin <- 0.002 # 0.002 
# gsMax	220
x1$paramsTranspiration$Gswmax <- 0.220  # 0.220

# Conductances
x1$paramsTranspiration$VCleaf_kmax <- 1/(1/k_SLApo + 1/k_LSym)
# x1$paramsTranspiration$VCleaf_kmax <- k_SLApo
rs_stem <- (1/x1$paramsTranspiration$VCstem_kmax)
rs_root <- (1/x1$paramsTranspiration$VCroot_kmax)
rs_stem_new <- (1/k_RSApo)*(rs_stem/(rs_stem+rs_root))
rs_root_new <- (1/k_RSApo)*(rs_root/(rs_stem+rs_root))
sum(rs_stem_new+rs_root_new)
x1$paramsTranspiration$VCstem_kmax <- 1/rs_stem_new
x1$paramsTranspiration$VCroot_kmax <- 1/rs_root_new
1/(1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)
x1$belowLayers$VCroot_kmax <- x1$belowLayers$VCroot_kmax*(rs_root/rs_root_new)
x1$paramsTranspiration$Plant_kmax <- 1/(1/x1$paramsTranspiration$VCleaf_kmax + 1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)

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

control$TPhase_gmin	<- 37.5
control$Q10_1_gmin <- 1.2
control$Q10_2_gmin <- 4.8

# gCrown0	150
control$gCrown0 <- 0.150
# gsMax	220
# gsNight	20
control$gs_NightFrac <- 20/220
# JarvisPAR	0.003
control$JarvisPAR  <- 0.003

# k_SSymInit	0.26
control$k_SSym <- 0.26
control$k_LSym <- k_LSym

# fFRBToLeaf	1
# fTRBToLeaf	0.8
control$fTRBToLeaf <- 0.8

# C_LApoInit	0.00001
# C_SApoInit	0.00002
control$C_LApoInit <- 0.00001
control$C_SApoInit <- 0.00002

#0.003
control$JarvisPAR <- 0.003

#Initialize input
x2 <- forest2spwbInput(pue_forest, pue_soil, SpParamsMED, control)
# Regulates the proportion of sunlit and shade leaves
x2$paramsInterception$kDIR <- 0.4
# K	0.5
x2$paramsInterception$kPAR <- 0.5
# canopyStorageParam	1.5
x2$paramsInterception$g <- 1.5/2.2
# gmin20	2
# gmin_S	2
x2$paramsTranspiration$Gswmin <- 0.002
# gsMax	220
x2$paramsTranspiration$Gswmax <- 0.300 ##0.220
# Tgs_sens	17
x2$paramsTranspiration$Gs_Tsens <- 17
# Tgs_optim	25
x2$paramsTranspiration$Gs_Toptim <- 25
# Name	Value
# Foliage	Evergreen
# P50_VC_Leaf	-6.4
x2$paramsTranspiration$VCleaf_P50 <- -6.4
# slope_VC_Leaf	30
x2$paramsTranspiration$VCleaf_slope <- 30
# P50_VC_Stem	-6.4
x2$paramsTranspiration$VCstem_P50 <- -6.4
x2$paramsTranspiration$VCroot_P50 <- -6.4
# slope_VC_Stem	30
x2$paramsTranspiration$VCstem_slope <- 30
x2$paramsTranspiration$VCroot_slope <- 30
P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
wb <- hydraulics_psi2Weibull(-6.4, P88)
x2$paramsTranspiration$VCleaf_c <- wb["c"]
x2$paramsTranspiration$VCleaf_d <- wb["d"]
x2$paramsTranspiration$VCstem_c <- x2$paramsTranspiration$VCleaf_c
x2$paramsTranspiration$VCroot_c <- x2$paramsTranspiration$VCleaf_c
x2$paramsTranspiration$VCstem_d <- x2$paramsTranspiration$VCleaf_d
x2$paramsTranspiration$VCroot_d <- x2$paramsTranspiration$VCleaf_d
# Conductances
x2$paramsTranspiration$VCleaf_kmax <- k_SLApo 
rs_stem <- (1/x2$paramsTranspiration$VCstem_kmax)
rs_root <- (1/x2$paramsTranspiration$VCroot_kmax)
rs_stem_new <- (1/k_RSApo)*(rs_stem/(rs_stem+rs_root))
rs_root_new <- (1/k_RSApo)*(rs_root/(rs_stem+rs_root))
sum(rs_stem_new+rs_root_new)
x2$paramsTranspiration$VCstem_kmax <- 1/rs_stem_new
x2$paramsTranspiration$VCroot_kmax <- 1/rs_root_new
1/(1/x2$paramsTranspiration$VCstem_kmax + 1/x2$paramsTranspiration$VCroot_kmax)
x2$belowLayers$VCroot_kmax <- x2$belowLayers$VCroot_kmax*(rs_root/rs_root_new)
x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCleaf_kmax + 
                                        1/x2$paramsTranspiration$VCstem_kmax + 
                                        1/x2$paramsTranspiration$VCroot_kmax+
                                        1/k_LSym)

# P12_gs	-1
# P88_gs	-2.7
x2$paramsTranspiration$Gs_slope <- (88.0 - 12.0)/(2.7 - 1);
x2$paramsTranspiration$Gs_P50 <- -1.0 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope/25)
# epsilonSym_Leaf	15
x2$paramsWaterStorage$LeafEPS <- 15
# PiFullTurgor_Leaf	-2.5
x2$paramsWaterStorage$LeafPI0 <- -2.5
# apoFrac_Leaf	0.4
x2$paramsWaterStorage$LeafAF <- 0.4

# PiFullTurgor_Stem	-2.5
x2$paramsWaterStorage$StemPI0 <- -2.5
# epsilonSym_Stem	15
x2$paramsWaterStorage$StemEPS <- 15
# apoFrac_Stem	0.4
# symFrac_Stem	0.2
x2$paramsWaterStorage$StemAF <- 0.4
# vol_Stem	15
# x2$paramsWaterStorage$Vsapwood <- 15


S2 <- spwb(x2, pue_meteo, latitude = pue_latitude, elevation = pue_elevation)
# shinyplot(S2)