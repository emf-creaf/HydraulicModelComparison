library(traits4models)
prades_soil <- function(){
  pr_soil_df <- read.table("Data/Prades/PRADES_soilData.txt", sep="\t", header = TRUE)
  # pr_soil_df$widths[4]<- 3000
  pr_soil <- soil(pr_soil_df)
  pr_soil$rfc <- c(45, 70, 80, 95)
  return(pr_soil)
}

prades_forest <- function(){
  pr_forest <- emptyforest()
  pr_forest$treeData <- read.table("Data/Prades/PRADES_treeData.txt", sep="\t", header = TRUE)
  pr_forest$treeData$Z50 <- 200
  pr_forest$treeData$Z95 <- 1000
  return(pr_forest)
}

prades_input <- function(control) {
  data(SpParamsES)
  SpParams <- modifySpParams(SpParamsES, read.table("Data/Prades/PRADES_customParams.txt", sep="\t", header = TRUE), subsetSpecies = TRUE)
  pr_forest <- prades_forest()
  pr_soil <- prades_soil()
  
  ps <- 1
  qi <- 2
  
  k_plant <- c(1.0, 0.8)

  rf_leaf <- c(0.40, 0.60)
  rf_leaf_apo <- c(0.15, 0.20)
  rf_leaf_symp <- c(0.25, 0.40)
  rf_stem <- c(0.30, 0.25)
  rf_root <- c(0.30, 0.15)
  
  if(control$transpirationMode == "Granier") {
    x0 <- forest2spwbInput(pr_forest, pr_soil, SpParams, control)
    
    # P.latifolia	Qilex	Phalepensis
    # P50_VC_Leaf	 -6.5	-6.4	 -4.79
    # slope_VC_Leaf	11	30	46
    # P50_VC_Stem	 -6.5	-6.4	 -4.79
    # slope_VC_Stem	11	30	46
    # P88 <- -6.5 + log((100.0/88.0)-1.0)*(25.0/11)
    # wb <- hydraulics_psi2Weibull(-6.5, P88)
    # x0$paramsTranspiration$VCleaf_c[pl] <- wb["c"]
    # x0$paramsTranspiration$VCleaf_d[pl] <- wb["d"]
    # x0$paramsTranspiration$VCstem_c[pl] <- wb["c"]
    # x0$paramsTranspiration$VCstem_d[pl] <- wb["d"]
    P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
    wb <- hydraulics_psi2Weibull(-6.4, P88)
    x0$paramsTranspiration$VCleaf_c[qi] <- wb["c"]
    x0$paramsTranspiration$VCleaf_d[qi] <- wb["d"]
    x0$paramsTranspiration$VCstem_c[qi] <- wb["c"]
    x0$paramsTranspiration$VCstem_d[qi] <- wb["d"]
    P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    wb <- hydraulics_psi2Weibull(-4.79, P88)
    x0$paramsTranspiration$VCleaf_c[ps] <- wb["c"]
    x0$paramsTranspiration$VCleaf_d[ps] <- wb["d"]
    x0$paramsTranspiration$VCstem_c[ps] <- wb["c"]
    x0$paramsTranspiration$VCstem_d[ps] <- wb["d"]
    return(x0)
  } else if(control$transpirationMode == "Sperry") {
    
    #Initialize input
    x1 <- forest2spwbInput(pr_forest, pr_soil, SpParams, control)
    # # Regulates the proportion of sunlit and shade leaves
    # x1$paramsInterception$kDIR <- 0.4
    # x1$paramsInterception$kPAR <- 0.5
    # 
    # # P.latifolia	Qilex	Phalepensis
    # # P50_VC_Leaf	 -6.5	-6.4	 -4.79
    # # slope_VC_Leaf	11	30	46
    # # P50_VC_Stem	 -6.5	-6.4	 -4.79
    # # slope_VC_Stem	11	30	46
    # # P88 <- -6.5 + log((100.0/88.0)-1.0)*(25.0/11)
    # # wb <- hydraulics_psi2Weibull(-6.5, P88)
    # # x1$paramsTranspiration$VCleaf_c[pl] <- wb["c"]
    # # x1$paramsTranspiration$VCleaf_d[pl] <- wb["d"]
    # # x1$paramsTranspiration$VCstem_c[pl] <- wb["c"]
    # # x1$paramsTranspiration$VCstem_d[pl] <- wb["d"]
    # # x1$paramsTranspiration$VCroot_c[pl] <- wb["c"]
    # # x1$paramsTranspiration$VCroot_d[pl] <- wb["d"]
    # P88 <- -6.4 + log((100.0/88.0)-1.0)*(25.0/30)
    # wb <- hydraulics_psi2Weibull(-6.4, P88)
    # x1$paramsTranspiration$VCleaf_P50[qi] <- -6.4
    # x1$paramsTranspiration$VCleaf_slope[qi] <- 30
    # x1$paramsTranspiration$VCleaf_c[qi] <- wb["c"]
    # x1$paramsTranspiration$VCleaf_d[qi] <- wb["d"]
    # x1$paramsTranspiration$VCstem_P50[qi] <- -6.4
    # x1$paramsTranspiration$VCstem_slope[qi] <- 30
    # x1$paramsTranspiration$VCstem_c[qi] <- wb["c"]
    # x1$paramsTranspiration$VCstem_d[qi] <- wb["d"]
    # x1$paramsTranspiration$VCroot_P50[qi] <- -6.4
    # x1$paramsTranspiration$VCroot_slope[qi] <- 30
    # x1$paramsTranspiration$VCroot_c[qi] <- wb["c"]
    # x1$paramsTranspiration$VCroot_d[qi] <- wb["d"]
    # P88 <- -4.79 + log((100.0/88.0)-1.0)*(25.0/46)
    # wb <- hydraulics_psi2Weibull(-4.79, P88)
    # x1$paramsTranspiration$VCleaf_P50[ps] <- -4.79
    # x1$paramsTranspiration$VCleaf_slope[ps] <- 46
    # x1$paramsTranspiration$VCleaf_c[ps] <- wb["c"]
    # x1$paramsTranspiration$VCleaf_d[ps] <- wb["d"]
    # x1$paramsTranspiration$VCstem_P50[ps] <- -4.79
    # x1$paramsTranspiration$VCstem_slope[ps] <- 46
    # x1$paramsTranspiration$VCstem_c[ps] <- wb["c"]
    # x1$paramsTranspiration$VCstem_d[ps] <- wb["d"]
    # x1$paramsTranspiration$VCroot_P50[ps] <- -4.79
    # x1$paramsTranspiration$VCroot_slope[ps] <- 46
    # x1$paramsTranspiration$VCroot_c[ps] <- wb["c"]
    # x1$paramsTranspiration$VCroot_d[ps] <- wb["d"]
    # 
    # # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # # epsilonSym_Leaf	12.38	12.38	15	5.31	5
    # # PiFullTurgor_Leaf	 -2.7	 -2.13	-2.5	 -1.5	-1.65
    # # apoFrac_Leaf	0.5	0.5	0.4	0.6	0.6
    # # x1$paramsWaterStorage$LeafEPS[pl] <- 12.38 
    # # x1$paramsWaterStorage$LeafPI0[pl] <- -2.13
    # # x1$paramsWaterStorage$LeafAF[pl] <- 0.5
    # x1$paramsWaterStorage$LeafEPS[qi] <- 15 
    # x1$paramsWaterStorage$LeafPI0[qi] <- -2.5
    # x1$paramsWaterStorage$LeafAF[qi] <- 0.4
    # x1$paramsWaterStorage$LeafEPS[ps] <- 5.31 
    # x1$paramsWaterStorage$LeafPI0[ps] <- -1.5
    # x1$paramsWaterStorage$LeafAF[ps] <- 0.6
    # 
    # # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # # PiFullTurgor_Stem	12.38	12.38	-2.5	-1.65	-1.65
    # # epsilonSym_Stem	 -2.7	 -2.13	15	5	5
    # # apoFrac_Stem	0.4	0.4	0.4	0.4	0.4
    # # symFrac_Stem	0.2	0.2	0.2	0.2	0.2
    # # x1$paramsWaterStorage$StemEPS[pl] <- 12.38 
    # # x1$paramsWaterStorage$StemPI0[pl] <- -2.13
    # # x1$paramsWaterStorage$StemAF[pl] <- 0.4
    # x1$paramsWaterStorage$StemEPS[qi] <- 15 
    # x1$paramsWaterStorage$StemPI0[qi] <- -2.5
    # x1$paramsWaterStorage$StemAF[qi] <- 0.4
    # x1$paramsWaterStorage$StemEPS[ps] <- 5
    # x1$paramsWaterStorage$StemPI0[ps] <- -1.65
    # x1$paramsWaterStorage$StemAF[ps] <- 0.4
    # 
    # Conductances
    x1$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*rf_leaf_symp)
    x1$paramsTranspiration$VCleafapo_kmax <- 1.0/((1.0/k_plant)*rf_leaf_apo)
    x1$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*rf_leaf)
    x1$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*rf_stem)
    x1$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*rf_root)
    for(i in c(qi, ps)) x1$belowLayers$VCroot_kmax[i,] <- x1$belowLayers$VCroot_kmax[i, ]*x1$paramsTranspiration$VCroot_kmax[i]/sum(x1$belowLayers$VCroot_kmax[i, ])
    x1$paramsTranspiration$Plant_kmax <- 1/(1/x1$paramsTranspiration$VCleaf_kmax + 1/x1$paramsTranspiration$VCstem_kmax + 1/x1$paramsTranspiration$VCroot_kmax)
    x1$paramsTranspiration$FR_leaf <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCleaf_kmax
    x1$paramsTranspiration$FR_stem <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCstem_kmax
    x1$paramsTranspiration$FR_root <- x1$paramsTranspiration$Plant_kmax/x1$paramsTranspiration$VCroot_kmax
    # 
    # # P. angustifolia	P.latifolia	Qilex	Phalepensis	Phalepensis_OLD
    # # gsMax	150	200	220	217.500963190747	217.500963190747
    # # gmin20	1.5	2	2	1	1.5
    # # x1$paramsTranspiration$Gswmax[pl] <- 0.220 
    # x1$paramsTranspiration$Gswmax[qi] <- 0.220 
    # x1$paramsTranspiration$Gswmax[ps] <- 0.2175 
    # # x1$paramsTranspiration$Gswmin[pl] <- 0.002 
    # x1$paramsTranspiration$Gswmin[qi] <- 0.002 
    # x1$paramsTranspiration$Gswmin[ps] <- 0.001 
    # 
    medfate:::.updateBelow(x1)
    return(x1)
  } else if(control$transpirationMode == "Sureau") {

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
    x2 <- forest2spwbInput(pr_forest, pr_soil, SpParams, control)
    # Regulates the proportion of sunlit and shade leaves
    x2$paramsInterception$kDIR <- 0.4
    x2$paramsInterception$kPAR <- 0.5


    #           Qilex	Psylvestris
    # P50_VC_Leaf	 -6.4	 -3.2
    # slope_VC_Leaf	30	17
    # P50_VC_Stem	 -6.4	 -3.2
    # slope_VC_Stem	30	17

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

    x2$paramsTranspiration$VCleaf_P50[ps] <- -3.2
    x2$paramsTranspiration$VCleaf_slope[ps] <- 17
    x2$paramsTranspiration$VCstem_P50[ps] <- -3.2
    x2$paramsTranspiration$VCstem_slope[ps] <- 17
    x2$paramsTranspiration$VCroot_P50[ps] <- -3.2
    x2$paramsTranspiration$VCroot_slope[ps] <- 17
    P88 <- -3.2 + log((100.0/88.0)-1.0)*(25.0/17)
    wb <- hydraulics_psi2Weibull(-3.2, P88)
    x2$paramsTranspiration$VCleaf_c[ps] <- wb["c"]
    x2$paramsTranspiration$VCleaf_d[ps] <- wb["d"]
    x2$paramsTranspiration$VCstem_c[ps] <- wb["c"]
    x2$paramsTranspiration$VCstem_d[ps] <- wb["d"]
    x2$paramsTranspiration$VCroot_c[ps] <- wb["c"]
    x2$paramsTranspiration$VCroot_d[ps] <- wb["d"]

    #                 Qilex	Psylvestris
    # epsilonSym_Leaf	  15	5.31
    # PiFullTurgor_Leaf	 -2.5	 -1.5
    # apoFrac_Leaf	    0.4	0.6	
    x2$paramsWaterStorage$LeafEPS[qi] <- 15
    x2$paramsWaterStorage$LeafPI0[qi] <- -2.5
    x2$paramsWaterStorage$LeafAF[qi] <- 0.4
    x2$paramsWaterStorage$LeafEPS[ps] <- 5.31
    x2$paramsWaterStorage$LeafPI0[ps] <- -1.5
    x2$paramsWaterStorage$LeafAF[ps] <- 0.6

    #                 Qilex	Psylvestris
    # PiFullTurgor_Stem	-2.5	-1.65
    # epsilonSym_Stem	  15	5
    # apoFrac_Stem	  0.4	0.4
    # symFrac_Stem	0.2	0.2
    x2$paramsWaterStorage$StemEPS[qi] <- 15
    x2$paramsWaterStorage$StemPI0[qi] <- -2.5
    x2$paramsWaterStorage$StemAF[qi] <- 0.4
    x2$paramsWaterStorage$StemEPS[ps] <- 5
    x2$paramsWaterStorage$StemPI0[ps] <- -1.65
    x2$paramsWaterStorage$StemAF[ps] <- 0.4

    # Conductances
    x2$paramsTranspiration$kleaf_symp <- 1.0/((1.0/k_plant)*rf_leaf_symp)
    x2$paramsTranspiration$VCleafapo_kmax <- 1.0/((1.0/k_plant)*rf_leaf_apo)
    x2$paramsTranspiration$VCleaf_kmax <- 1.0/((1.0/k_plant)*rf_leaf)
    x2$paramsTranspiration$VCstem_kmax <- 1.0/((1.0/k_plant)*rf_stem)
    x2$paramsTranspiration$VCroot_kmax <- 1.0/((1.0/k_plant)*rf_root)
    for(i in c(qi, ps)) x2$belowLayers$VCroot_kmax[i,] <- x2$belowLayers$VCroot_kmax[i, ]*x2$paramsTranspiration$VCroot_kmax[i]/sum(x2$belowLayers$VCroot_kmax[i, ])
    x2$paramsTranspiration$Plant_kmax <- 1/(1/x2$paramsTranspiration$VCleaf_kmax +
                                              1/x2$paramsTranspiration$VCstem_kmax +
                                              1/x2$paramsTranspiration$VCroot_kmax)
    x2$paramsTranspiration$FR_leaf <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCleaf_kmax
    x2$paramsTranspiration$FR_stem <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCstem_kmax
    x2$paramsTranspiration$FR_root <- x2$paramsTranspiration$Plant_kmax/x2$paramsTranspiration$VCroot_kmax

    #       Qilex  P. sylvestris
    # gsMax	  220	180
    # gmin20	2  1
    x2$paramsTranspiration$Gswmax[qi] <- 0.220
    x2$paramsTranspiration$Gswmax[ps] <- 0.180
    x2$paramsTranspiration$Gswmin[qi] <- 0.002
    x2$paramsTranspiration$Gswmin[ps] <- 0.001



    #         Qilex	Psylvestris
    # P12_gs	-1      -1.36
    # P88_gs	-2.7   -2.14
    x2$paramsTranspiration$Gs_slope[qi] <- (88.0 - 12.0)/(2.7 - 1);
    x2$paramsTranspiration$Gs_P50[qi] <- -1.0 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope[qi]/25)
    x2$paramsTranspiration$Gs_slope[ps] <- (88.0 - 12.0)/(2.14 - 1.36);
    x2$paramsTranspiration$Gs_P50[ps] <- -1.36 + log(0.12/0.88)/(x2$paramsTranspiration$Gs_slope[ps]/25)

    medfate:::.updateBelow(x2)
    return(x2)
  }
}
