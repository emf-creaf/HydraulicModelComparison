# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)
library(tidyverse)
library(cowplot)

source("Rscripts/Ancillary.R")

# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270
pue_slope <- 0
pue_aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo[1:200,]
meteo$DOY <- 200
meteo$JulianDay <- 200
meteo$MinTemperature <- 20
meteo$MaxTemperature <- 30
meteo$MinRelativeHumidity <- 30
meteo$MaxRelativeHumidity <- 90
meteo$WindSpeed <- 2
meteo$Radiation  <- 30
meteo$Precipitation <- 0 # To force dessication

# Sperry simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$sapFluidityVariation <- FALSE
control$leafCavitationEffects <- TRUE
control$rhizosphereOverlap <- "total"

#Initialize input
x1 <- puechabon_input(control)

#Change canopy and soil variables
x1$canopy$Tair <- 29
x1$canopy$Cair <- 386
x1$canopy$VPair <- 1.688
x1$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1 <- spwb(x1, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)


# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$sapFluidityVariation <- FALSE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE
control$rhizosphereOverlap <- "total"


x2 <- puechabon_input(control)

#Change canopy and soil variables
x2$canopy$Tair <- 29
x2$canopy$Cair <- 386
x2$canopy$VPair <- 1.688
x2$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S2 <- spwb(x2, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)



# Sperry (segmented) ------------------------------------------------------
x3 <- x1
wb <- hydraulics_psi2Weibull(x2$paramsTranspiration$Gs_P50, -2.7)
x3$paramsTranspiration$VCleaf_c <- wb["c"]
x3$paramsTranspiration$VCleaf_d <- wb["d"]
x3$control$leafCavitationEffects <- TRUE
S3 <- spwb(x3, meteo, 
           latitude = pue_latitude, elevation = pue_elevation, 
           slope = pue_slope, aspect = pue_aspect)


# Plots -------------------------------------------------------------------

res1 <- extractSubdaily(S1, output = "LeafPsi") |>
  dplyr::mutate(Model = "Sperry (non-segmented)") |>
  dplyr::rename(LeafPsi = T1_168)
res2 <- extractSubdaily(S2, output = "LeafPsi") |>
  dplyr::mutate(Model = "Sureau") |>
  dplyr::rename(LeafPsi = T1_168)
res3 <- extractSubdaily(S3, output = "LeafPsi") |>
  dplyr::mutate(Model = "Sperry (leaf segmented)") |>
  dplyr::rename(LeafPsi = T1_168)
res <- dplyr::bind_rows(res1, res2, res3) |>
  dplyr::filter(!is.na(LeafPsi))
ggplot(res, aes(x=datetime, y = LeafPsi))+
  geom_line(aes(col = Model))+
  xlab("")+ ylab("Leaf water potential (MPa)")+
  theme_bw()

extractDailyVar <- function(var) {
  res1 <- as.data.frame(S1$Plants[[var]])
  res1$date <- as.Date(row.names(res1))
  res1$Model <- "Sperry (non-segmented)"
  row.names(res1) <- NULL
  names(res1)[1] <- "var"
  res2 <- as.data.frame(S2$Plants[[var]])
  res2$date <- as.Date(row.names(res2))
  res2$Model <- "Sureau"
  row.names(res2) <- NULL
  names(res2)[1] <- "var"
  res3 <- as.data.frame(S3$Plants[[var]])
  res3$date <- as.Date(row.names(res3))
  res3$Model <- "Sperry (leaf segmented)"
  row.names(res3) <- NULL
  names(res3)[1] <- "var"
  res <- dplyr::bind_rows(res1, res2, res3) |>
    dplyr::filter(!is.na(var))
  return(res)
}

res <- extractDailyVar("Transpiration")
ggplot(res, aes(x=date, y = var))+
  geom_line(aes(col = Model))+
  xlab("")+ ylab("Transpiration (mm)")+
  theme_bw()

res <- extractDailyVar("LeafPLC")
res$var <- res$var*100
ggplot(res, aes(x=date, y = var))+
  geom_line(aes(col = Model))+
  xlab("")+ ylab("Leaf PLC (%)")+
  theme_bw()

res <- extractDailyVar("StemPLC")
res$var <- res$var*100
ggplot(res, aes(x=date, y = var))+
  geom_line(aes(col = Model))+
  xlab("")+ ylab("Stem PLC (%)")+
  theme_bw()

res <- extractDailyVar("dEdP")
res$var <- res$var
ggplot(res, aes(x=date, y = var))+
  geom_line(aes(col = Model))+
  xlab("")+ ylab("k_plant")+
  theme_bw()


# Plots -------------------------------------------------------------------
p1 <- plot(S1, "SoilPsi")+ylim(c(-5,0))+labs(title="Sperry")
p2 <- plot(S2, "SoilPsi")+ylim(c(-5,0))+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/SoilPsi_Puechabon_Dessication.png", p, width = 6, height = 8)

p1 <- plot(S1, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "LeafPsiRange", bySpecies = TRUE)+ylim(c(-5,0))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/LeafPsiRange_Puechabon_Dessication.png", p, width = 6, height = 8)

p1 <- plot(S1, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "Transpiration", bySpecies = TRUE)+ylim(c(0,3.5))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/Transpiration_Puechabon_Dessication.png", p, width = 6, height = 8)


p1 <- plot(S1, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "GSWMax_SL", bySpecies = TRUE)+ylim(c(0,0.3))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/StomatalConductance_Puechabon_Dessication.png", p, width = 6, height = 8)

p1 <- plot(S1, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "StemPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/StemPLC_Puechabon_Dessication.png", p, width = 6, height = 8)

p1 <- plot(S1, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "LeafPLC", bySpecies = TRUE)+ylim(c(0,100))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/LeafPLC_Puechabon_Dessication.png", p, width = 6, height = 8)

p1 <- plot(S1, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sperry")
p2 <- plot(S2, "SoilPlantConductance", bySpecies = TRUE)+ylim(c(0,1))+theme(legend.position = "none")+labs(title="Sureau")
p <-plot_grid(p1, p2,
              nrow = 2)
ggsave2("Plots/Puechabon_Dessication/SoilPlantConductance_Puechabon_Dessication.png", p, width = 6, height = 8)

