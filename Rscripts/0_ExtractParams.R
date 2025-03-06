
extract_spparams<-function(x, param_selection = c( "Name","Z50", "Z95", "SLA",
                                                   "Al2As", "r635",  
                                                   "Vmax298","Jmax298", 
                                                   "maxFMC","LeafPI0", "LeafEPS", "LeafAF", "StemAF",
                                                   "Gs_P50", "Gswmin", 
                                                   "Plant_kmax", "VCleaf_kmax", "VCleafapo_kmax", "kleaf_symp", "VCstem_kmax", "VCroot_kmax",
                                                   "VCstem_c", "VCstem_d","VCleaf_c", "VCleaf_d","VCroot_c", "VCroot_d",
                                                   "VCroot_P50", "VCroot_slope","VCstem_P50", "VCstem_slope","VCleaf_P50", "VCleaf_slope",
                                                   "Kmax_stemxylem")) {
  df <-cbind(as.data.frame(x$cohorts),
             as.data.frame(x$below),
             as.data.frame(x$paramsAnatomy),
             as.data.frame(x$paramsTranspiration),
             as.data.frame(x$paramsWaterStorage))
  df <- df |>
    dplyr::select(any_of(param_selection))
  row.names(df) <- NULL
  return(df)
}
