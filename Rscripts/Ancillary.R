puechabon_soil<- function() {
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
  return(pue_soil)
}

fontblanche_soil <- function(){
  fb_soil_df <- defaultSoilParams(n = 3)
  fb_soil_df$widths <- c(200, 800, 3000)
  fb_soil_df$rfc <- c(50, 70, 94)
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

yatir_soil<- function() {
  pue_soil_df <- defaultSoilParams(n = 3)
  pue_soil_df$widths <- c(200, 300, 3000)
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
  return(pue_soil)
}

puechabon_forest <- function(){
  pue_forest <- emptyforest(1)
  pue_forest$treeData$Species <- "Quercus ilex"
  pue_forest$treeData$DBH <- 9.11
  pue_forest$treeData$Height <- 530.2
  pue_forest$treeData$N <- 1750
  pue_forest$treeData$Z50 <- 200
  pue_forest$treeData$Z95 <- 1000
  pue_forest$treeData$LAI <- 2.2
  return(pue_forest)
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

fontblanche_forest <- function(){
  fb_forest <- emptyforest()
  fb_forest$treeData <- read.table("Data/FontBlanche/FONBLA_treeData.txt", sep="\t", header = TRUE)
  fb_forest$treeData$Z50 <- 200
  fb_forest$treeData$Z95 <- 1000
  fb_forest$treeData$LAI <- species_LAI(fb_forest, SpParamsMED)
  # remove Phillyrea
  # fb_forest$treeData <- fb_forest$treeData[-1, ]
  fb_forest$treeData$LAI <- fb_forest$treeData$LAI*(2.7/sum(fb_forest$treeData$LAI))
  return(fb_forest)
}

puechabon_input <- function(control) {
  data(SpParamsMED)
  pue_forest <- puechabon_forest()
  pue_soil <- puechabon_soil()
  
  k_plant = 0.8
  k_RSApo1 <- 0.364
  k_RSApo2 <- 0.903
  k_RSApo3 <- 0.732
  # k_RSApo <- k_RSApo1 + k_RSApo2 + k_RSApo3
  k_RSApo <- 2
  k_SLApo <- 4
  k_LSym <- 2
  1/(1/k_RSApo + 1/k_SLApo + 1/k_LSym)
  
  P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
  wb <- hydraulics_psi2Weibull(-6.4, P88)
  
  if(control$transpirationMode == "Granier") {
    x0 <- forest2spwbInput(pue_forest, pue_soil, SpParamsMED, control)
    x0$paramsTranspiration$VCstem_c <- wb["c"]
    x0$paramsTranspiration$VCstem_d <- wb["d"]
    return(x0)
  }
  else if(control$transpirationMode == "Sperry") {
    #Initialize input
    x1 <- forest2spwbInput(pue_forest, pue_soil, SpParamsMED, control)
    
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
    # medfate:::.updateBelow(x1)
    return(x1)
  }
  else if(control$transpirationMode =="Cochard"){
    control$TPhase_gmin	<- 37.5
    control$Q10_1_gmin <- 1.2
    control$Q10_2_gmin <- 4.8
    
    # gCrown0	150
    control$gCrown0 <- 0.150
    # gsMax	220
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
    x2$paramsTranspiration$Gswmax <- 0.220 ##0.220
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
    # medfate:::.updateBelow(x2)
    return(x2)
  }
  return(NULL)
}

fontblanche_input <- function(control) {
  data(SpParamsMED)
  fb_forest <- fontblanche_forest()
  fb_soil <- fontblanche_soil()
  
  pl <- 1
  ph <- 2
  qi <- 3
  
  k_plant <- rep(0, 3)
  k_plant[pl] <- 1.55
  k_plant[qi] <- 0.8
  k_plant[ph] <- 0.55
  
  if(control$transpirationMode == "Granier") {
    x0 <- forest2spwbInput(fb_forest, fb_soil, SpParamsMED, control)
    
    # P.latifolia	Qilex	Phalepensis
    # P50_VC_Leaf	 -6.5	-6.4	 -4.79
    # slope_VC_Leaf	11	30	46
    # P50_VC_Stem	 -6.5	-6.4	 -4.79
    # slope_VC_Stem	11	30	46
    P88 <- -6.5 + log((100.0/88.0)-1.0)*(25.0/11)
    wb <- hydraulics_psi2Weibull(-6.5, P88)
    x0$paramsTranspiration$VCleaf_c[pl] <- wb["c"]
    x0$paramsTranspiration$VCleaf_d[pl] <- wb["d"]
    x0$paramsTranspiration$VCstem_c[pl] <- wb["c"]
    x0$paramsTranspiration$VCstem_d[pl] <- wb["d"]
    P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
    wb <- hydraulics_psi2Weibull(-6.4, P88)
    x0$paramsTranspiration$VCleaf_c[qi] <- wb["c"]
    x0$paramsTranspiration$VCleaf_d[qi] <- wb["d"]
    x0$paramsTranspiration$VCstem_c[qi] <- wb["c"]
    x0$paramsTranspiration$VCstem_d[qi] <- wb["d"]
    P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    wb <- hydraulics_psi2Weibull(-4.79, P88)
    x0$paramsTranspiration$VCleaf_c[ph] <- wb["c"]
    x0$paramsTranspiration$VCleaf_d[ph] <- wb["d"]
    x0$paramsTranspiration$VCstem_c[ph] <- wb["c"]
    x0$paramsTranspiration$VCstem_d[ph] <- wb["d"]
    return(x0)
  }
  else if(control$transpirationMode == "Sperry") {
    
    #Initialize input
    x1 <- forest2spwbInput(fb_forest, fb_soil, SpParamsMED, control)
    # Regulates the proportion of sunlit and shade leaves
    x1$paramsInterception$kDIR <- 0.4
    x1$paramsInterception$kPAR <- 0.5
    
    # P.latifolia	Qilex	Phalepensis
    # P50_VC_Leaf	 -6.5	-6.4	 -4.79
    # slope_VC_Leaf	11	30	46
    # P50_VC_Stem	 -6.5	-6.4	 -4.79
    # slope_VC_Stem	11	30	46
    P88 <- -6.5 + log((100.0/88.0)-1.0)*(25.0/11)
    wb <- hydraulics_psi2Weibull(-6.5, P88)
    x1$paramsTranspiration$VCleaf_c[pl] <- wb["c"]
    x1$paramsTranspiration$VCleaf_d[pl] <- wb["d"]
    x1$paramsTranspiration$VCstem_c[pl] <- wb["c"]
    x1$paramsTranspiration$VCstem_d[pl] <- wb["d"]
    x1$paramsTranspiration$VCroot_c[pl] <- wb["c"]
    x1$paramsTranspiration$VCroot_d[pl] <- wb["d"]
    P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
    wb <- hydraulics_psi2Weibull(-6.4, P88)
    x1$paramsTranspiration$VCleaf_c[qi] <- wb["c"]
    x1$paramsTranspiration$VCleaf_d[qi] <- wb["d"]
    x1$paramsTranspiration$VCstem_c[qi] <- wb["c"]
    x1$paramsTranspiration$VCstem_d[qi] <- wb["d"]
    x1$paramsTranspiration$VCroot_c[qi] <- wb["c"]
    x1$paramsTranspiration$VCroot_d[qi] <- wb["d"]
    P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    wb <- hydraulics_psi2Weibull(-4.79, P88)
    x1$paramsTranspiration$VCleaf_c[ph] <- wb["c"]
    x1$paramsTranspiration$VCleaf_d[ph] <- wb["d"]
    x1$paramsTranspiration$VCstem_c[ph] <- wb["c"]
    x1$paramsTranspiration$VCstem_d[ph] <- wb["d"]
    x1$paramsTranspiration$VCroot_c[ph] <- wb["c"]
    x1$paramsTranspiration$VCroot_d[ph] <- wb["d"]
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # epsilonSym_Leaf	12.38	12.38	15	5.31	5
    # PiFullTurgor_Leaf	 -2.7	 -2.13	-2.5	 -1.5	-1.65
    # apoFrac_Leaf	0.5	0.5	0.4	0.6	0.6
    x1$paramsWaterStorage$LeafEPS[pl] <- 12.38 
    x1$paramsWaterStorage$LeafPI0[pl] <- -2.13
    x1$paramsWaterStorage$LeafAF[pl] <- 0.5
    x1$paramsWaterStorage$LeafEPS[qi] <- 15 
    x1$paramsWaterStorage$LeafPI0[qi] <- -2.5
    x1$paramsWaterStorage$LeafAF[qi] <- 0.4
    x1$paramsWaterStorage$LeafEPS[ph] <- 5.31 
    x1$paramsWaterStorage$LeafPI0[ph] <- -1.5
    x1$paramsWaterStorage$LeafAF[ph] <- 0.6
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # PiFullTurgor_Stem	12.38	12.38	-2.5	-1.65	-1.65
    # epsilonSym_Stem	 -2.7	 -2.13	15	5	5
    # apoFrac_Stem	0.4	0.4	0.4	0.4	0.4
    # symFrac_Stem	0.2	0.2	0.2	0.2	0.2
    x1$paramsWaterStorage$StemEPS[pl] <- 12.38 
    x1$paramsWaterStorage$StemPI0[pl] <- -2.13
    x1$paramsWaterStorage$StemAF[pl] <- 0.4
    x1$paramsWaterStorage$StemEPS[qi] <- 15 
    x1$paramsWaterStorage$StemPI0[qi] <- -2.5
    x1$paramsWaterStorage$StemAF[qi] <- 0.4
    x1$paramsWaterStorage$StemEPS[ph] <- 5
    x1$paramsWaterStorage$StemPI0[ph] <- -1.65
    x1$paramsWaterStorage$StemAF[ph] <- 0.4
    
    # Conductances
    # resistances 40% leaf / 30% stem / 30% root
    x1$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*0.4)
    x1$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*0.3)
    x1$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*0.3)
    for(i in c(pl, qi, ph)) x1$belowLayers$VCroot_kmax[i,] <- x1$belowLayers$VCroot_kmax[i, ]*x1$paramsTranspiration$VCroot_kmax[i]/sum(x1$belowLayers$VCroot_kmax[i, ])
    x1$paramsTranspiration$Plant_kmax <- 1/(1/x1$paramsTranspiration$VCleaf_kmax + 1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # gsMax	150	200	220	217.500963190747	217.500963190747
    # gmin20	1.5	2	2	1	1.5
    x1$paramsTranspiration$Gswmax[pl] <- 0.220 
    x1$paramsTranspiration$Gswmax[qi] <- 0.220 
    x1$paramsTranspiration$Gswmax[ph] <- 0.2175 
    x1$paramsTranspiration$Gswmin[pl] <- 0.002 
    x1$paramsTranspiration$Gswmin[qi] <- 0.002 
    x1$paramsTranspiration$Gswmin[ph] <- 0.001 
    
    medfate:::.updateBelow(x1)
    return(x1)
  }
  else if(control$transpirationMode == "Cochard") {
    
    control$TPhase_gmin	<- 37.5
    control$Q10_1_gmin <- 1.2
    control$Q10_2_gmin <- 4.8
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
    x2 <- forest2spwbInput(fb_forest, fb_soil, SpParamsMED, control)
    # Regulates the proportion of sunlit and shade leaves
    x2$paramsInterception$kDIR <- 0.4
    x2$paramsInterception$kPAR <- 0.5
    
    
    # P.latifolia	Qilex	Phalepensis
    # P50_VC_Leaf	 -6.5	-6.4	 -4.79
    # slope_VC_Leaf	11	30	46
    # P50_VC_Stem	 -6.5	-6.4	 -4.79
    # slope_VC_Stem	11	30	46
    x2$paramsTranspiration$VCleaf_P50[pl] <- -6.5
    x2$paramsTranspiration$VCleaf_slope[pl] <- 11
    x2$paramsTranspiration$VCstem_P50[pl] <- -6.5
    x2$paramsTranspiration$VCstem_slope[pl] <- 11
    x2$paramsTranspiration$VCroot_P50[pl] <- -6.5
    x2$paramsTranspiration$VCroot_slope[pl] <- 11
    P88 <- -6.5 + log((100.0/88.0)-1.0)*(25.0/11)
    wb <- hydraulics_psi2Weibull(-6.5, P88)
    x2$paramsTranspiration$VCleaf_c[pl] <- wb["c"]
    x2$paramsTranspiration$VCleaf_d[pl] <- wb["d"]
    x2$paramsTranspiration$VCstem_c[pl] <- wb["c"]
    x2$paramsTranspiration$VCstem_d[pl] <- wb["d"]
    x2$paramsTranspiration$VCroot_c[pl] <- wb["c"]
    x2$paramsTranspiration$VCroot_d[pl] <- wb["d"]
    
    x2$paramsTranspiration$VCleaf_P50[qi] <- -6.4
    x2$paramsTranspiration$VCleaf_slope[qi] <- 30
    x2$paramsTranspiration$VCstem_P50[qi] <- -6.4
    x2$paramsTranspiration$VCstem_slope[qi] <- 30
    x2$paramsTranspiration$VCroot_P50[qi] <- -6.4
    x2$paramsTranspiration$VCroot_slope[qi] <- 30
    P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
    wb <- hydraulics_psi2Weibull(-6.4, P88)
    x2$paramsTranspiration$VCleaf_c[qi] <- wb["c"]
    x2$paramsTranspiration$VCleaf_d[qi] <- wb["d"]
    x2$paramsTranspiration$VCstem_c[qi] <- wb["c"]
    x2$paramsTranspiration$VCstem_d[qi] <- wb["d"]
    x2$paramsTranspiration$VCroot_c[qi] <- wb["c"]
    x2$paramsTranspiration$VCroot_d[qi] <- wb["d"]
    
    x2$paramsTranspiration$VCleaf_P50[ph] <- -4.79
    x2$paramsTranspiration$VCleaf_slope[ph] <- 46
    x2$paramsTranspiration$VCstem_P50[ph] <- -4.79
    x2$paramsTranspiration$VCstem_slope[ph] <- 46
    x2$paramsTranspiration$VCroot_P50[ph] <- -4.79
    x2$paramsTranspiration$VCroot_slope[ph] <- 46
    P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    wb <- hydraulics_psi2Weibull(-4.79, P88)
    x2$paramsTranspiration$VCleaf_c[ph] <- wb["c"]
    x2$paramsTranspiration$VCleaf_d[ph] <- wb["d"]
    x2$paramsTranspiration$VCstem_c[ph] <- wb["c"]
    x2$paramsTranspiration$VCstem_d[ph] <- wb["d"]
    x2$paramsTranspiration$VCroot_c[ph] <- wb["c"]
    x2$paramsTranspiration$VCroot_d[ph] <- wb["d"]
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # epsilonSym_Leaf	12.38	12.38	15	5.31	5
    # PiFullTurgor_Leaf	 -2.7	 -2.13	-2.5	 -1.5	-1.65
    # apoFrac_Leaf	0.5	0.5	0.4	0.6	0.6
    x2$paramsWaterStorage$LeafEPS[pl] <- 12.38 
    x2$paramsWaterStorage$LeafPI0[pl] <- -2.13
    x2$paramsWaterStorage$LeafAF[pl] <- 0.5
    x2$paramsWaterStorage$LeafEPS[qi] <- 15 
    x2$paramsWaterStorage$LeafPI0[qi] <- -2.5
    x2$paramsWaterStorage$LeafAF[qi] <- 0.4
    x2$paramsWaterStorage$LeafEPS[ph] <- 5.31 
    x2$paramsWaterStorage$LeafPI0[ph] <- -1.5
    x2$paramsWaterStorage$LeafAF[ph] <- 0.6
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # PiFullTurgor_Stem	12.38	12.38	-2.5	-1.65	-1.65
    # epsilonSym_Stem	 -2.7	 -2.13	15	5	5
    # apoFrac_Stem	0.4	0.4	0.4	0.4	0.4
    # symFrac_Stem	0.2	0.2	0.2	0.2	0.2
    x2$paramsWaterStorage$StemEPS[pl] <- 12.38 
    x2$paramsWaterStorage$StemPI0[pl] <- -2.13
    x2$paramsWaterStorage$StemAF[pl] <- 0.4
    x2$paramsWaterStorage$StemEPS[qi] <- 15 
    x2$paramsWaterStorage$StemPI0[qi] <- -2.5
    x2$paramsWaterStorage$StemAF[qi] <- 0.4
    x2$paramsWaterStorage$StemEPS[ph] <- 5
    x2$paramsWaterStorage$StemPI0[ph] <- -1.65
    x2$paramsWaterStorage$StemAF[ph] <- 0.4
    
    # Conductances
    # resistances 20% leaf symp / 20% leaf apo / 30% stem / 30% root
    x2$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*0.20)
    x2$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*0.20)
    x2$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*0.3)
    x2$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*0.3)
    for(i in c(pl, qi, ph)) x2$belowLayers$VCroot_kmax[i,] <- x2$belowLayers$VCroot_kmax[i, ]*x2$paramsTranspiration$VCroot_kmax[i]/sum(x2$belowLayers$VCroot_kmax[i, ])
    x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCleaf_kmax + 
                                              1/x2$paramsTranspiration$VCstem_kmax + 
                                              1/x2$paramsTranspiration$VCroot_kmax+
                                              1/x2$paramsTranspiration$kleaf_symp)
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # gsMax	150	200	220	217.500963190747	217.500963190747
    # gmin20	1.5	2	2	1	1.5
    x2$paramsTranspiration$Gswmax[pl] <- 0.220 
    x2$paramsTranspiration$Gswmax[qi] <- 0.220 
    x2$paramsTranspiration$Gswmax[ph] <- 0.2175 
    x2$paramsTranspiration$Gswmin[pl] <- 0.002 
    x2$paramsTranspiration$Gswmin[qi] <- 0.002 
    x2$paramsTranspiration$Gswmin[ph] <- 0.001 
    
    
    
    # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # P12_gs	-2.2	 -1.65	-1	-1.36	-1.36
    # P88_gs	-3.8	 -2.5	-2.7	 -2.14	-2.33
    x2$paramsTranspiration$Gs_slope[pl] <- (88.0 - 12.0)/(2.5 - 1.65);
    x2$paramsTranspiration$Gs_P50[pl] <- -1.65 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope[pl]/25)
    x2$paramsTranspiration$Gs_slope[qi] <- (88.0 - 12.0)/(2.7 - 1);
    x2$paramsTranspiration$Gs_P50[qi] <- -1.0 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope[qi]/25)
    x2$paramsTranspiration$Gs_slope[ph] <- (88.0 - 12.0)/(2.14 - 1.36);
    x2$paramsTranspiration$Gs_P50[ph] <- -1.36 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope[ph]/25)
    
    medfate:::.updateBelow(x2)
    return(x2)
  }
}

yatir_input <- function(control) {
  data(SpParamsMED)
  yat_forest <- yatir_forest()
  yat_soil <- yatir_soil()
  
  k_plant <- 0.55
  k_RSApo <- 2
  k_SLApo <- 4
  k_LSym <- 2
  1/(1/k_RSApo + 1/k_SLApo + 1/k_LSym)
  
  P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
  wb <- hydraulics_psi2Weibull(-4.79, P88)
  
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
    x1$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*0.4)
    x1$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*0.3)
    x1$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*0.3)
    x1$belowLayers$VCroot_kmax <- x1$belowLayers$VCroot_kmax*x1$paramsTranspiration$VCroot_kmax/sum(x1$belowLayers$VCroot_kmax)
    x1$paramsTranspiration$Plant_kmax <- 1/(1/x1$paramsTranspiration$VCleaf_kmax + 1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)
    
    medfate:::.updateBelow(x1)
    return(x1)
  }
  else if(control$transpirationMode =="Cochard"){
    control$TPhase_gmin	<- 37.5
    control$Q10_1_gmin <- 1.2
    control$Q10_2_gmin <- 4.8
    
    # gCrown0	150
    control$gCrown0 <- 0.150
    # gsMax	220
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
    # P50_VC_Leaf	-4.79
    x2$paramsTranspiration$VCleaf_P50 <- -4.79
    # slope_VC_Leaf	46
    x2$paramsTranspiration$VCleaf_slope <- 46
    # P50_VC_Stem	-4.79
    x2$paramsTranspiration$VCstem_P50 <- -4.79
    x2$paramsTranspiration$VCroot_P50 <- -4.79
    # slope_VC_Stem	46
    x2$paramsTranspiration$VCstem_slope <- 46
    x2$paramsTranspiration$VCroot_slope <- 46
    
    P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    wb <- hydraulics_psi2Weibull(-4.79, P88)
    
    x2$paramsTranspiration$VCleaf_c <- wb["c"]
    x2$paramsTranspiration$VCleaf_d <- wb["d"]
    x2$paramsTranspiration$VCstem_c <- x2$paramsTranspiration$VCleaf_c
    x2$paramsTranspiration$VCroot_c <- x2$paramsTranspiration$VCleaf_c
    x2$paramsTranspiration$VCstem_d <- x2$paramsTranspiration$VCleaf_d
    x2$paramsTranspiration$VCroot_d <- x2$paramsTranspiration$VCleaf_d
    # Conductances
    # resistances 20% leaf symp / 20% leaf apo / 30% stem / 30% root
    x2$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*0.20)
    x2$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*0.20)
    x2$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*0.3)
    x2$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*0.3)
    x2$belowLayers$VCroot_kmax <- x2$belowLayers$VCroot_kmax*x2$paramsTranspiration$VCroot_kmax/sum(x2$belowLayers$VCroot_kmax)
    x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCleaf_kmax + 
                                              1/x2$paramsTranspiration$VCstem_kmax + 
                                              1/x2$paramsTranspiration$VCroot_kmax+
                                              1/x2$paramsTranspiration$kleaf_symp)
    
    # P12_gs	-1
    # P88_gs	-2.7
    x1$paramsTranspiration$Gs_slope <- (88.0 - 12.0)/(2.14 - 1.36);
    x1$paramsTranspiration$Gs_P50 <- -1.36 + log(0.12/0.88)/(x1$paramsTranspiration$Gs_slope/25)
    
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
