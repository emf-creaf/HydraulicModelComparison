yatir_soil<- function() {
  yat_soil_df <- defaultSoilParams(n = 3)
  yat_soil_df$widths <- c(200, 300, 3000)
  yat_soil_df$rfc <- c(75, 82, 94)
  yat_soil_df$clay <- rep(39, 3)
  yat_soil_df$sand <- rep(26, 3)
  yat_soil_df$om <- c(6,3,1)
  yat_soil_df$bd <-rep(1.45, 3)
  
  yat_soil <- soil(yat_soil_df)
  yat_soil$VG_theta_sat <- rep(0.5, 3)
  yat_soil$VG_theta_res <- rep(0.098, 3)
  yat_soil$VG_n <- rep(1.55, 3)
  # 1 cm = 0.00009804139432 MPa
  yat_soil$VG_alpha <- rep(0.0005/0.00009804139432, 3)
  yat_soil$Ksat <- rep(655.2934*1.69, 3)
  return(yat_soil)
}

yatir_forest <- function(){
  yat_forest <- emptyforest(1)
  yat_forest$treeData$Species <- "Pinus halepensis"
  yat_forest$treeData$DBH <- 17
  yat_forest$treeData$Height <- 1000
  yat_forest$treeData$N <- 300
  yat_forest$treeData$Z50 <- 200
  yat_forest$treeData$Z95 <- 1000
  yat_forest$treeData$LAI <- 1.5
  return(yat_forest)
}

yatir_input <- function(control) {
  data(SpParamsMED)
  yat_forest <- yatir_forest()
  yat_soil <- yatir_soil()
  
  k_plant <- 0.55
  
  rf_leaf <- 0.40
  rf_leaf_apo <- 0.15
  rf_leaf_symp <- 0.25
  rf_stem <- 0.30
  rf_root <- 0.30

  P50 <- -4.79
  slope <- 46
  P88 <- P50 + log((100.0/88.0)-1.0)*(25.0/slope)
  wb <- hydraulics_psi2Weibull(P50, P88)
  
  if(control$transpirationMode == "Granier") {
    x0 <- forest2spwbInput(yat_forest, yat_soil, SpParamsMED, control)
    x0$paramsTranspiration$VCstem_c <- wb["c"]
    x0$paramsTranspiration$VCstem_d <- wb["d"]
    return(x0)
  }
  else if(control$transpirationMode == "Sperry") {
    #Initialize input
    x1 <- forest2spwbInput(yat_forest, yat_soil, SpParamsMED, control)
    
    x1$paramsInterception$kDIR <- 0.4
    # K	0.5
    x1$paramsInterception$kPAR <- 0.5
    # canopyStorageParam	1.5
    x1$paramsInterception$g <- 1.5/2.2
    # P50_VC_Leaf	-6.4
    # slope_VC_Leaf	30
    # wbl <- hydraulics_psi2Weibull(-4.5, -5.5)
    # wbl <- hydraulics_psi2Weibull(-3.5, -4.5)
    wbl <- wb
    x1$paramsTranspiration$VCleaf_c <- wbl["c"]
    x1$paramsTranspiration$VCstem_c <- wb["c"]
    x1$paramsTranspiration$VCroot_c <- wb["c"]
    x1$paramsTranspiration$VCleaf_d <- wbl["d"]
    x1$paramsTranspiration$VCstem_d <- wb["d"]
    x1$paramsTranspiration$VCroot_d <- wb["d"]
    
    # gmin20	2
    # gmin_S	2
    x1$paramsTranspiration$Gswmin <- 0.001 # 0.002 
    # gsMax	220
    x1$paramsTranspiration$Gswmax <- 0.2175  # 0.2175
    
    # Conductances
    x1$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*rf_leaf_symp)
    x1$paramsTranspiration$VCleafapo_kmax <- 1.0/((1.0/k_plant)*rf_leaf_apo)
    x1$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*rf_leaf)
    x1$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*rf_stem)
    x1$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*rf_root)
    x1$belowLayers$VCroot_kmax <- x1$belowLayers$VCroot_kmax*x1$paramsTranspiration$VCroot_kmax/sum(x1$belowLayers$VCroot_kmax)
    x1$paramsTranspiration$Plant_kmax <- 1/(1/x1$paramsTranspiration$VCleaf_kmax + 1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)
    x1$paramsTranspiration$FR_leaf <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCleaf_kmax
    x1$paramsTranspiration$FR_stem <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCstem_kmax
    x1$paramsTranspiration$FR_root <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCroot_kmax
    
    medfate:::.updateBelow(x1)
    return(x1)
  }
  else if(control$transpirationMode =="Sureau"){

    # gCrown0	150
    control$gCrown0 <- 0.150
    # gsMax	220
    # JarvisPAR	0.003
    control$JarvisPAR  <- 0.003
    
    # k_SSymInit	0.26
    control$k_SSym <- 0.26

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
    x2 <- forest2spwbInput(yat_forest, yat_soil, SpParamsMED, control)
    # Regulates the proportion of sunlit and shade leaves
    x2$paramsInterception$kDIR <- 0.4
    # K	0.5
    x2$paramsInterception$kPAR <- 0.5
    # canopyStorageParam	1.5
    x2$paramsInterception$g <- 1.5/2.2
    # gmin20	2
    # gmin_S	2
    x2$paramsTranspiration$Gswmin <- 0.001
    # gsMax	220
    x2$paramsTranspiration$Gswmax <- 0.2175 ##0.2175
    # Tgs_sens	17
    x2$paramsTranspiration$Gs_Tsens <- 17
    # Tgs_optim	25
    x2$paramsTranspiration$Gs_Toptim <- 25
    # Name	Value
    # Foliage	Evergreen
    x2$paramsTranspiration$VCleaf_P50 <- P50
    x2$paramsTranspiration$VCstem_P50 <- P50
    x2$paramsTranspiration$VCroot_P50 <- P50
    x2$paramsTranspiration$VCleaf_slope <- slope
    x2$paramsTranspiration$VCstem_slope <- slope
    x2$paramsTranspiration$VCroot_slope <- slope
    
    x2$paramsTranspiration$VCleaf_c <- wb["c"]
    x2$paramsTranspiration$VCleaf_d <- wb["d"]
    x2$paramsTranspiration$VCstem_c <- x2$paramsTranspiration$VCleaf_c
    x2$paramsTranspiration$VCroot_c <- x2$paramsTranspiration$VCleaf_c
    x2$paramsTranspiration$VCstem_d <- x2$paramsTranspiration$VCleaf_d
    x2$paramsTranspiration$VCroot_d <- x2$paramsTranspiration$VCleaf_d
    # Conductances
    # resistances 20% leaf symp / 20% leaf apo / 30% stem / 30% root
    x2$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*rf_leaf_symp)
    x2$paramsTranspiration$VCleafapo_kmax <- 1.0/((1.0/k_plant)*rf_leaf_apo)
    x2$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*rf_leaf)
    x2$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*rf_stem)
    x2$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*rf_root)
    x2$belowLayers$VCroot_kmax <- x2$belowLayers$VCroot_kmax*x2$paramsTranspiration$VCroot_kmax/sum(x2$belowLayers$VCroot_kmax)
    x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCleaf_kmax + 
                                              1/x2$paramsTranspiration$VCstem_kmax + 
                                              1/x2$paramsTranspiration$VCroot_kmax)
    x2$paramsTranspiration$FR_leaf <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCleaf_kmax
    x2$paramsTranspiration$FR_stem <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCstem_kmax
    x2$paramsTranspiration$FR_root <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCroot_kmax
    
    # P12_gs	-1.47
    # P88_gs	-2.17
    x2$paramsTranspiration$Gs_slope <- (88.0 - 12.0)/(2.17 - 1.47);
    x2$paramsTranspiration$Gs_P50 <- -1.47 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope/25)
    
    # epsilonSym_Leaf	5.3
    x2$paramsWaterStorage$LeafEPS <- 5.31
    # PiFullTurgor_Leaf	-1.5
    x2$paramsWaterStorage$LeafPI0 <- -1.5
    # apoFrac_Leaf	0.6
    x2$paramsWaterStorage$LeafAF <- 0.6
    
    # PiFullTurgor_Stem	-1.5
    x2$paramsWaterStorage$StemPI0 <- -1.5
    # epsilonSym_Stem	5.3
    x2$paramsWaterStorage$StemEPS <- 5.31
    # apoFrac_Stem	0.6
    x2$paramsWaterStorage$StemAF <- 0.6
    medfate:::.updateBelow(x2)
    return(x2)
  }
  return(NULL)
}
