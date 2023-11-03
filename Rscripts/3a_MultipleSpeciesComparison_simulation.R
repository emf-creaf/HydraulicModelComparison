library(readr)
library(medfate)
data(SpParamsMED)

# K option 1
# Use k_plant as final value and fixed resistance fractions (40 leaves, 30 stem, 40 roots)
# K option 2
# Use k_plant to estimate k_leaf and max_height to estimate k_stem, then roots account for 40%


# Terrain -----------------------------------------------------------------
pue_latitude <- 43.74139
pue_elevation <- 270
pue_slope <- 0
pue_aspect <- 0

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


common_soil <- function(){
  fb_soil_df <- defaultSoilParams(n = 3)
  fb_soil_df$widths <- c(200, 800, 3000)
  fb_soil_df$rfc <- c(30, 50, 94.5)
  fb_soil_df$clay <- rep(39, 3)
  fb_soil_df$sand <- rep(26, 3)
  fb_soil_df$om <- c(6,3,1)
  fb_soil_df$bd <-rep(1.45, 3)
  
  fb_soil <- soil(fb_soil_df)
  fb_soil$VG_theta_sat <- rep(0.5, 3)
  fb_soil$VG_theta_res <- rep(0.098, 3)
  fb_soil$VG_n <- rep(1.55, 3)
  # 1 cm = 0.00009804139432 MPa
  fb_soil$VG_alpha <- rep(0.0005/0.00009804139432, 3)
  fb_soil$Ksat <- rep(655.2934*1.69, 3)
  
  return(fb_soil)
}

create_forest <- function(sp_name, height_cm){
  pue_forest <- emptyforest(1)
  pue_forest$treeData$Species <- sp_name
  pue_forest$treeData$DBH <- 25
  pue_forest$treeData$Height <- height_cm
  pue_forest$treeData$N <- NA
  pue_forest$treeData$Z50 <- 306.5
  pue_forest$treeData$Z95 <- 1000
  pue_forest$treeData$LAI <- 2.0
  pue_forest$treeData$CR <- 0.5
  return(pue_forest)
}


# Species parameters ------------------------------------------------------

NewParams <- read_delim("Data/NewParamsOK_First.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE) |> as.data.frame()


# K_leaf from Wolfe et al.
NewParams$fRleaf <- pmax(0.2, 0.70 - 0.011*NewParams$HeightMax)
NewParams$kleaf <- 1/(NewParams$fRleaf/NewParams$kplantConsensus)

# Gs_P50 and Gs_slope
NewParams$Gs_slope <- (88-12)/(abs(NewParams$Pgs88) - abs(NewParams$Pgs12))
# NewParams$Gs_P50 <- NewParams$Pgs88 - (25/NewParams$Gs_slope)*log((1/0.88)-1)
NewParams$Gs_P50 <-(NewParams$Pgs88 + NewParams$Pgs12)/2.0
1/(1+exp((NewParams$Gs_slope[1]/25)*(NewParams$Pgs12[1]- (NewParams$Gs_P50[1]))))
1/(1+exp((NewParams$Gs_slope[1]/25)*(NewParams$Pgs88[1]- (NewParams$Gs_P50[1]))))

# P12_VC and P88_VC
NewParams$P12_VC <- NewParams$P50_VC + log((1/0.12)-1.0)*(25.0/NewParams$Slope_VC)
NewParams$P88_VC <- NewParams$P50_VC + log((1/0.88)-1.0)*(25.0/NewParams$Slope_VC)

# Nleaf
NewParams$Nleaf <- 13



SpParams <- medfateutils::initSpParams(NewParams$Species, SpParamsDefinition)
SpParams$GrowthForm <- ifelse(NewParams$GrowthForm=="T", "Tree", "Shrub")
SpParams$LifeForm <- "Phanerophyte"
SpParams$LeafShape <- "Broad"
SpParams$LeafSize <- "Large"
SpParams$PhenologyType <- "continuous-evergreen"
SpParams$Gswmax <- NewParams$gs_consensus
SpParams$Gswmin <- NewParams$gmin_mmol/1000
SpParams$Kmax_stemxylem <- 1 # NewParams$Kstem[sp_index]
SpParams$Hmax <- NewParams$HeightMax*100
SpParams$Hmed <- SpParams$Hmax*0.7
SpParams$SLA <- 1000/NewParams$LMA_gm2
SpParams$Al2As <- 1000/NewParams$HV
SpParams$VCstem_P12 <- NewParams$P12_VC
SpParams$VCstem_P88 <- NewParams$P88_VC
SpParams$VCstem_P50 <- NewParams$P50_VC
SpParams$VCstem_slope <- NewParams$Slope_VC
SpParams$VCleaf_P12 <- NewParams$P12_VC
SpParams$VCleaf_P88 <- NewParams$P88_VC
SpParams$VCleaf_P50 <- NewParams$P50_VC
SpParams$VCleaf_slope <- NewParams$Slope_VC
SpParams$VCleaf_kmax <- NewParams$kleaf
SpParams$VCroot_P12 <- NewParams$P12_VC
SpParams$VCroot_P88 <- NewParams$P88_VC
SpParams$VCroot_P50 <- NewParams$P50_VC
SpParams$VCroot_slope <- NewParams$Slope_VC
SpParams$LeafPI0 <- NewParams$Pi0
SpParams$LeafEPS <- NewParams$Esymp
SpParams$LeafAF <- NA
SpParams$StemPI0 <- NA
SpParams$StemEPS <- NA
SpParams$StemAF <- NA
SpParams$Gs_P50 <- NewParams$Gs_P50
SpParams$Gs_slope <- NewParams$Gs_slope
SpParams$Nleaf <- NewParams$Nleaf
SpParams$Vmax298 <- NA
SpParams$Jmax298 <- NA
write.table(SpParams, "Tables/SpParams.txt", sep="\t", quote=FALSE)



# Simulations -------------------------------------------------------------
soil <- common_soil()

nspp <- nrow(SpParams)

results <- data.frame(Species = SpParams$Name)
results$actual_kplant <- NA
results$psi_min <- NA
results$E_max <- NA
results$day_closure <- NA
results$day_failure_leaf <- NA
results$day_failure_stem <- NA

simSperry <- TRUE
simSperry_segmented <- TRUE
simSureau <- TRUE

results_sperry <- results
results_sperry_segmented <-results
results_sureau <- results

for(sp_index in 1:nspp) {
  cat(paste0(sp_index, " Species name: ", results$Species[sp_index]))

  cat(paste0(" FOREST "))
  forest <- create_forest( NewParams$Species[sp_index],  NewParams$HeightMax[sp_index]*100)
  

  if(simSperry) {
    tryCatch({
      cat(paste0(" SPERRY "))
      #Initialize control parameters
      control <- defaultControl("Sperry")
      control$subdailyResults <- TRUE
      control$cavitationRefill <- "none"
      control$bareSoilEvaporation <- FALSE
      control$sapFluidityVariation <- FALSE
      control$leafCavitationEffects <- TRUE
      control$rhizosphereOverlap <- "total"
      control$verbose <- FALSE
      
      #Initialize input
      x1 <- forest2spwbInput(forest, soil, SpParams, control)
      
      #Change canopy and soil variables
      x1$canopy$Tair <- 29
      x1$canopy$Cair <- 386
      x1$canopy$VPair <- 1.688
      x1$soil$Temp <- c(32,29,27.71661)
      #Call simulation function
      S1 <- spwb(x1, meteo,
                 latitude = pue_latitude, elevation = pue_elevation,
                 slope = pue_slope, aspect = pue_aspect)
      
      results_sperry$actual_kplant[sp_index] <- S1$spwbOutput$paramsTranspiration$Plant_kmax[1]
      results_sperry$psi_min[sp_index] <- S1$Plants$LeafPsiMin[1]
      results_sperry$E_max[sp_index] <- max(S1$subdaily[[1]]$SunlitLeavesInst$E)
      
      gswmax <- S1$SunlitLeaves$GSWMax
      gswmax <- gswmax[!is.na(gswmax)]
      
      results_sperry$day_closure[sp_index] <- which(gswmax < x1$paramsTranspiration$Gswmax*0.05)[1] # First day with maximum gs < 5% of gsmax
      
      results_sperry$day_failure_leaf[sp_index] <- which(is.na(S1$Plants$LeafPLC))[1]
      results_sperry$day_failure_stem[sp_index] <- which(is.na(S1$Plants$StemPLC))[1]
      
      if(is.na(results_sperry$day_closure[sp_index])) {
        results_sperry$day_closure[sp_index] = results_sperry$day_failure_leaf[sp_index]
      }
    })
  }
  
  if(simSperry_segmented) {
    tryCatch({
      cat(paste0(" SPERRY SEGMENTED "))
      #Initialize control parameters
      control <- defaultControl("Sperry")
      control$subdailyResults <- TRUE
      control$cavitationRefill <- "none"
      control$bareSoilEvaporation <- FALSE
      control$sapFluidityVariation <- FALSE
      control$leafCavitationEffects <- TRUE
      control$rhizosphereOverlap <- "total"
      control$verbose <- FALSE

      #Initialize input
      SpParamsOne_S <- SpParams[sp_index,]
      SpParamsOne_S$VCleaf_P50 <- SpParamsOne_S$Gs_P50
      SpParamsOne_S$VCleaf_slope <- SpParamsOne_S$Gs_slope
      SpParamsOne_S$VCleaf_P12 <- NA
      SpParamsOne_S$VCleaf_P88 <- NA
      x1 <- forest2spwbInput(forest, soil, SpParamsOne_S, control)

      #Change canopy and soil variables
      x1$canopy$Tair <- 29
      x1$canopy$Cair <- 386
      x1$canopy$VPair <- 1.688
      x1$soil$Temp <- c(32,29,27.71661)
      #Call simulation function
      S1 <- spwb(x1, meteo,
                 latitude = pue_latitude, elevation = pue_elevation,
                 slope = pue_slope, aspect = pue_aspect)
      
      results_sperry_segmented$actual_kplant[sp_index] <- S1$spwbOutput$paramsTranspiration$Plant_kmax[1]
      results_sperry_segmented$psi_min[sp_index] <- S1$Plants$LeafPsiMin[1]
      results_sperry_segmented$E_max[sp_index] <- max(S1$subdaily[[1]]$SunlitLeavesInst$E)
      
      gswmax <- S1$SunlitLeaves$GSWMax
      gswmax <- gswmax[!is.na(gswmax)]
      
      results_sperry_segmented$day_closure[sp_index] <- which(gswmax < x1$paramsTranspiration$Gswmax*0.05)[1] # First day with maximum gs < 5% of gsmax
      
      results_sperry_segmented$day_failure_leaf[sp_index] <- which(is.na(S1$Plants$LeafPLC))[1]
      results_sperry_segmented$day_failure_stem[sp_index] <- which(is.na(S1$Plants$StemPLC))[1]
      if(is.na(results_sperry_segmented$day_closure[sp_index])) {
        results_sperry_segmented$day_closure[sp_index] = results_sperry_segmented$day_failure_leaf[sp_index]
      }
    })
  }
  
  if(simSureau) {
    #########################################################  
    # SUREAU SIMULATIONS
    # No cavitation refill
    # Leaf cuticular transpiration (for consistency with Sperry's use of gmin) but not stem cuticular transpiration
    # Plant capacitance disabled
    # No bare soil evaporation
    # Nocturnal transpiration equal to leaf cuticular transpiration (i.e. Gs_night = 0)
    # Leaf hydraulic conductance divided between apoplasmic (10% resistance) and symplasmic (90% resistance)
    #########################################################  
    tryCatch({
      cat(paste0(" SUREAU "))
      control <- defaultControl("Cochard")
      control$subdailyResults <- TRUE
      control$cavitationRefill <- "none"
      control$bareSoilEvaporation <- FALSE
      control$plantCapacitance <- TRUE
      control$sapFluidityVariation <- FALSE
      control$leafCuticularTranspiration <- TRUE
      control$stemCuticularTranspiration <- FALSE
      control$rhizosphereOverlap <- "total"
      control$stomatalSubmodel <- "Baldocchi"
      control$gs_NightFrac <- 0.01
      control$verbose <- FALSE
      
      x2 <- forest2spwbInput(forest, soil, SpParams, control)
      kleaf_tot <- x2$paramsTranspiration$VCleaf_kmax
      x2$paramsTranspiration$VCleaf_kmax <- 1/(0.1*(1/kleaf_tot))
      x2$paramsTranspiration$kleaf_symp <- 1/(0.9*(1/kleaf_tot))
      x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCroot_kmax +1/x2$paramsTranspiration$VCstem_kmax +1/x2$paramsTranspiration$VCleaf_kmax + 1/x2$paramsTranspiration$kleaf_symp)
      
      #Change canopy and soil variables
      x2$canopy$Tair <- 29
      x2$canopy$Cair <- 386
      x2$canopy$VPair <- 1.688
      x2$soil$Temp <- c(32,29,27.71661)
      #Call simulation function
      S2 <- spwb(x2, meteo, 
                 latitude = pue_latitude, elevation = pue_elevation, 
                 slope = pue_slope, aspect = pue_aspect)
      
      results_sureau$actual_kplant[sp_index] <- S2$spwbOutput$paramsTranspiration$Plant_kmax[1]
      results_sureau$psi_min[sp_index] <- S2$Plants$LeafPsiMin[1]
      results_sureau$E_max[sp_index] <- max(S2$subdaily[[1]]$SunlitLeavesInst$E)
      
      gswmax <- S2$SunlitLeaves$GSWMax
      gswmax <- gswmax[!is.na(gswmax)]
      
      results_sureau$day_closure[sp_index] <- which(gswmax < x2$paramsTranspiration$Gswmax*0.05)[1] # First day with maximum gs < 5% of gsmax
      
      results_sureau$day_failure_leaf[sp_index] <- which(is.na(S2$Plants$LeafPLC))[1]
      results_sureau$day_failure_stem[sp_index] <- which(is.na(S2$Plants$StemPLC))[1]
      if(is.na(results_sureau$day_closure[sp_index])) {
        results_sureau$day_closure[sp_index] = results_sureau$day_failure_leaf[sp_index]
      }
    })
  }
  
  cat(paste0("\n"))
}

if(simSperry) write.table(results_sperry, "Tables/results_sperry.txt", sep="\t", quote=FALSE)
if(simSperry_segmented) write.table(results_sperry_segmented, "Tables/results_sperry_segmented.txt", sep="\t", quote=FALSE)
if(simSureau) write.table(results_sureau, "Tables/results_sureau.txt", sep="\t", quote=FALSE)


