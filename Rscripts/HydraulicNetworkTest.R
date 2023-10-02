# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(ggplot2)


# Forest/topo preparation ------------------------------------------------------
# Example with one single cohort of Q. ilex
data("exampleforestMED2")
forest <- exampleforestMED2
forest$treeData <- forest$treeData[2, ,drop = FALSE]
forest$treeData$LAI <- 2.5
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
meteo$MinTemperature <- meteo$MinTemperature+10
meteo$MaxTemperature <- meteo$MaxTemperature+10
meteo$Precipitation <- 0 # To force dessication


# Soil definition ---------------------------------------------------------
df_soil <- defaultSoilParams(3)
df_soil$widths <- c(300,1200,2500)
df_soil$rfc <- c(40,65,97)
soil <- soil(df_soil)

stableSoilW <- matrix(1, nrow = nrow(meteo), ncol = nrow(df_soil))

# Sperry simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Sperry")
control$subdailyResults <- TRUE
# control$cavitationRefill <- "rate"
control$cavitationRefill <- "none"

#Initialize input
x1 <- forest2spwbInput(forest, soil, SpParamsMED, control)

#Call simulation function
S1 <- spwb(x1, meteo, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)
S1s <- pwb(x1, meteo, W = stableSoilW, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)


# Sureau simulation -------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
#Initialize input
x2 <- forest2spwbInput(forest, soil, SpParamsMED, control)

cn <- initCochardNetworks(x2)[[1]]
#Call simulation function
S2 <- spwb(x2, meteo, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)
S2s <- pwb(x2, meteo, W = stableSoilW, 
           latitude = latitude, elevation = elevation, 
           slope = slope, aspect = aspect)

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
plot(S1s, "SoilPlantConductance")+theme(legend.position = "none")
plot(S2s, "SoilPlantConductance")+theme(legend.position = "none")

plot(S1, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))
plot(S2, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))
plot(S1s, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))
plot(S2s, "StemPsi", subdaily = TRUE)+theme(legend.position = "none")+ylim(c(-8,0))

plot(S1, "StemPLC") +theme(legend.position = "none")
plot(S2, "StemPLC") +theme(legend.position = "none")
plot(S1s, "StemPLC") +theme(legend.position = "none")
plot(S2s, "StemPLC") +theme(legend.position = "none")

plot(S1, "LeafPLC") +theme(legend.position = "none")
plot(S2, "LeafPLC") +theme(legend.position = "none")
plot(S1s, "LeafPLC") +theme(legend.position = "none")
plot(S2s, "LeafPLC") +theme(legend.position = "none")

plot(S1, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "StemRWC", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafRWC") +theme(legend.position = "none")
plot(S2, "LeafRWC") +theme(legend.position = "none")
plot(S1s, "LeafRWC") +theme(legend.position = "none")
plot(S2s, "LeafRWC") +theme(legend.position = "none")

plot(S1, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))
plot(S2, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))
plot(S1s, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))
plot(S2s, "LeafPsiRange")+theme(legend.position = "none")+ylim(c(-5,0))


plot(S1, "SoilPsi") +theme(legend.position = "none")
plot(S2, "SoilPsi") +theme(legend.position = "none")
plot(S1s, "SoilPsi") +theme(legend.position = "none")
plot(S2s, "SoilPsi") +theme(legend.position = "none")

# Photosynthesis drops with stomatal closure in SUREAU
plot(S1, "PlantGrossPhotosynthesis")+theme(legend.position = "none")
plot(S2, "PlantGrossPhotosynthesis")+theme(legend.position = "none")
plot(S1s, "PlantGrossPhotosynthesis")+theme(legend.position = "none")
plot(S2s, "PlantGrossPhotosynthesis")+theme(legend.position = "none")

plot(S1, "PlantTranspiration")+theme(legend.position = "none")
plot(S2, "PlantTranspiration")+theme(legend.position = "none")
plot(S1s, "PlantTranspiration")+theme(legend.position = "none")
plot(S2s, "PlantTranspiration")+theme(legend.position = "none")


plot(S1, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafTemperature", subdaily = TRUE)+theme(legend.position = "none")

# Sureau differentiates more strongly between sunlit and shade leaves than Sperry (cause -> optimization)
plot(S1, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")

plot(S1$subdaily[[2]], "LeafStomatalConductance")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafStomatalConductance")+theme(legend.position = "none")
plot(S1$subdaily[[2]], "LeafTranspiration")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafTranspiration")+theme(legend.position = "none")
plot(S1$subdaily[[2]], "LeafPsi")+theme(legend.position = "none")
plot(S2$subdaily[[2]], "LeafPsi")+theme(legend.position = "none")

plot(S1, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafStomatalConductance", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafTranspiration", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafGrossPhotosynthesis", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafVPD", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafPsi", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafAbsorbedSWR", subdaily = TRUE)+theme(legend.position = "none")

plot(S1, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")
plot(S1s, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")
plot(S2s, "LeafAbsorbedPAR", subdaily = TRUE)+theme(legend.position = "none")



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
