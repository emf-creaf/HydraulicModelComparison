# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)


# Forest/topo preparation ------------------------------------------------------
data("exampleforestMED2")
forest <- exampleforestMED2
# forest$treeData <- forest$treeData[1, ,drop = FALSE] ## P. halepensis
forest$treeData <- forest$treeData[2, ,drop = FALSE] ## Q. ilex
# forest$treeData$LAI <- 2
forest$shrubData <- forest$shrubData[numeric(0), ,drop = FALSE]
forest$herbHeight <- NA
forest$herbLAI <- NA

latitude <- 41.82592
elevation <- 100
slope <- 0
aspect <- 0

# Weather preparation -----------------------------------------------------
data("examplemeteo")
meteo <- examplemeteo
meteo$DOY <- 200
meteo$JulianDay <- 200
meteo$MinTemperature <- 20
meteo$MaxTemperature <- 30
meteo$MinRelativeHumidity <- 30
meteo$MaxRelativeHumidity <- 90
meteo$WindSpeed <- 2
meteo$Radiation  <- 30
meteo$Precipitation <- 0 # To force dessication


# Soil definition ---------------------------------------------------------
df_soil <- defaultSoilParams(3)
df_soil$widths <- c(300,1200,2500)
df_soil$rfc <- c(40,65,97)
soil <- soil(df_soil)

stableSoilW <- matrix(1, nrow = nrow(meteo), ncol = nrow(df_soil))

SpParamsMED$VCleaf_P12 <- SpParamsMED$VCstem_P12
SpParamsMED$VCleaf_P50 <- SpParamsMED$VCstem_P50
SpParamsMED$VCleaf_P88 <- SpParamsMED$VCstem_P88
SpParamsMED$VCleaf_slope <- SpParamsMED$VCstem_slope

SpParamsMED$VCroot_P12 <- SpParamsMED$VCstem_P12
SpParamsMED$VCroot_P50 <- SpParamsMED$VCstem_P50
SpParamsMED$VCroot_P88 <- SpParamsMED$VCstem_P88
SpParamsMED$VCroot_slope <- SpParamsMED$VCstem_slope

# Sperry simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE


#Initialize input
x1 <- forest2spwbInput(forest, soil, SpParamsMED, control)
x1$paramsTranspiration$VCstem_kmax <- 2
x1$paramsTranspiration$VCleaf_kmax <- 0.7
x1$canopy$Tair <- 29
x1$Cair <- 386
x1$VPair <- 1.688
x1$soil$Temp <- c(32,29,27.71661)

#Call simulation function
S1 <- spwb(x1, meteo, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)




# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
control$bareSoilEvaporation <- FALSE
control$plantCapacitance <- TRUE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE

control$gs_NightFrac = 0 
#Initialize input
x2 <- forest2spwbInput(forest, soil, SpParamsMED, control)
x2$paramsTranspiration$Gs_P50 <- -2.5 # c(-2.5, -3.5)
x2$paramsTranspiration$Gs_slope <- 80
x2$paramsTranspiration$VCstem_kmax <- 2
x2$paramsTranspiration$VCleaf_kmax <- 1.5
x2$canopy$Tair <- 29
x2$Cair <- 386
x2$VPair <- 1.688
x2$soil$Temp <- c(32,29,27.71661)
x2$paramsTranspiration$Gswmin[1] <- 0.001
# 1/(1/0.7+1/2+1/1.5)

# cn <- initCochardNetworks(x2)

#Call simulation function
S2 <- spwb(x2, meteo, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)

# x3 <- x2
# x3$control$rootDisconnection <- TRUE
# S3 <- spwb(x3, meteo, 
#            latitude = latitude, elevation = elevation, 
#            slope = slope, aspect = aspect)

plot(S1$subdaily[[2]], "SoilPlantConductance")+theme(legend.position = "none")+ ylim(c(0,1))
plot(S2$subdaily[[2]], "SoilPlantConductance")+theme(legend.position = "none")+ ylim(c(0,1))
plot(S1$subdaily[[1]], "LeafTranspiration")+theme(legend.position = "none")
plot(S2$subdaily[[1]], "LeafTranspiration")+theme(legend.position = "none")
plot(S1$subdaily[[2]], "LeafGrossPhotosynthesis")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafGrossPhotosynthesis")+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafAbsorbedSWR")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafAbsorbedSWR")+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafStomatalConductance")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafStomatalConductance")+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafTemp")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafTemp")+theme(legend.position = "none")
plot(S1$subdaily[[2]], "LeafVPD")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafVPD")+theme(legend.position = "none")


plot(S1$subdaily[[2]], "LeafPLC")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafPLC")+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafRWC")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafRWC")+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafSympPsi")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafSympPsi")+theme(legend.position = "none")

plot(S1$subdaily[[1]], "StemPsi")+theme(legend.position = "none")
plot(S2$subdaily[[1]], "StemPsi")+theme(legend.position = "none")
plot(S1$subdaily[[1]], "LeafPsi")+theme(legend.position = "none")
plot(S2$subdaily[[1]], "LeafPsi")+theme(legend.position = "none")

plot(S2$subdaily[[1]], "StemRWC")+theme(legend.position = "none")

# Plots -------------------------------------------------------------------
hydraulics_vulnerabilityCurvePlot(x1, type = "leaf", relative = TRUE)+theme(legend.position = "none")
hydraulics_vulnerabilityCurvePlot(x2, type = "leaf", vulnerabilityFunction = "Sigmoid", relative = TRUE)+theme(legend.position = "none")
hydraulics_vulnerabilityCurvePlot(x1, type = "stem", relative = TRUE)+theme(legend.position = "none")
hydraulics_vulnerabilityCurvePlot(x2, type = "stem", vulnerabilityFunction = "Sigmoid", relative = TRUE)+theme(legend.position = "none")
hydraulics_vulnerabilityCurvePlot(x1, type = "root", relative = TRUE)+theme(legend.position = "none")
hydraulics_vulnerabilityCurvePlot(x2, type = "root", vulnerabilityFunction = "Sigmoid", relative = TRUE)+theme(legend.position = "none")


# Differences between both models
plot(S1, "SoilPlantConductance")+theme(legend.position = "none")
plot(S2, "SoilPlantConductance")+theme(legend.position = "none")
plot(S3, "SoilPlantConductance")+theme(legend.position = "none")

plot(S1, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))
plot(S2, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))

plot(S1, "StemPLC") +theme(legend.position = "none")
plot(S2, "StemPLC") +theme(legend.position = "none")
plot(S3, "StemPLC") +theme(legend.position = "none")

plot(S1, "LeafPLC") +theme(legend.position = "none")
plot(S2, "LeafPLC") +theme(legend.position = "none")
plot(S3, "LeafPLC") +theme(legend.position = "none")

plot(S1, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafRWC") +theme(legend.position = "none")
plot(S2, "LeafRWC") +theme(legend.position = "none")

plot(S1, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))
plot(S2, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))


plot(S1, "SoilPsi") +theme(legend.position = "none")
plot(S2, "SoilPsi") +theme(legend.position = "none")

plot(S1, "PlantGrossPhotosynthesis")+theme(legend.position = "none")
plot(S2, "PlantGrossPhotosynthesis")+theme(legend.position = "none")

plot(S1, "PlantTranspiration")+theme(legend.position = "none")+ylim(c(0,2.5))
plot(S2, "PlantTranspiration")+theme(legend.position = "none")+ylim(c(0,2.5))

plot(S1, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")


plot(S1, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")



SoilW <- S2$Soil[,1:3]
SoilPsi <- S2$Soil[,16:18]

df_list <- vector("list", nrow(SoilW))
for(day in 1:nrow(SoilW)) {
  date <- rownames(SoilW)[day]
  W <- as.numeric(SoilW[day, ])
  Einst<- as.numeric(S2$subdaily[[day]]$PlantsInst$E)
  LeafPsiInstSperry <- as.numeric(S1$subdaily[[day]]$PlantsInst$LeafPsi)
  StemPsiInstSperry <- as.numeric(S1$subdaily[[day]]$PlantsInst$StemPsi)
  RootPsiInst <- as.numeric(S1$subdaily[[day]]$PlantsInst$RootPsi)
  LeafPsiInstCochard <- as.numeric(S2$subdaily[[day]]$PlantsInst$LeafPsi)
  StemPLCInstCochard <- as.numeric(S2$subdaily[[day]]$PlantsInst$StemPLC)
  StemPsiInstCochard<- as.numeric(S2$subdaily[[day]]$PlantsInst$StemPsi)

  control <- defaultControl("Sperry")
  soil$W <- W
  xi <- forest2spwbInput(forest,soil, SpParamsMED, control)
  sn <- hydraulics_initSperryNetworks(xi)[[1]]
  # Soil psi is equal to previous day final soil water potential

  StemPsiInstSperry2 <- rep(NA, length(Einst))
  LeafPsiInstSperry2 <- rep(NA, length(Einst))
  for(i in 1:length(Einst)) {
    sn$PLCstem <- StemPLCInstCochard[i]
    # sol <- hydraulics_E2psiNetwork(Einst[i], sn)
    sol <- hydraulics_E2psiAboveground(Einst[i], RootPsiInst[i], sn)
    LeafPsiInstSperry2[i] <- sol$psiLeaf
    StemPsiInstSperry2[i] <- sol$psiStem
  }
  df_list[[day]] <- data.frame(date = rep(date, length(Einst)), step = 1:length(Einst),
                               RootPsiCochard = RootPsiInst,
                               StemPsiCochard = StemPsiInstCochard, 
                               StemPsiSperry1 = StemPsiInstSperry,
                               StemPsiSperry2 = StemPsiInstSperry2,
                               LeafPsiCochard = LeafPsiInstCochard,
                               LeafPsiSperry1 = LeafPsiInstSperry,
                               LeafPsiSperry2 = LeafPsiInstSperry2)
}
df <- dplyr::bind_rows(df_list)
df$date <- as.Date(df$date)
df$dateTime <- as.POSIXct(df$date) + 3600*(df$step-1)

ggplot(df)+
  geom_line(aes(x = dateTime, y = LeafPsiCochard), col = "red", alpha= 0.5)+
  geom_line(aes(x = dateTime, y = LeafPsiSperry1), col = "black", alpha= 0.5)+
  ylim(c(-5,0))


ggplot(df)+
  geom_line(aes(x = dateTime, y = StemPsiCochard), col = "red", alpha= 0.5)+
  geom_line(aes(x = dateTime, y = StemPsiSperry1), col = "black", alpha= 0.5)+
  ylim(c(-5,0))


