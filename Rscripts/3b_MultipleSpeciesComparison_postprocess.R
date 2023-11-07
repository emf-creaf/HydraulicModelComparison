library(medfate)
library(ggplot2)

# Process results ---------------------------------------------------------
SpParams <- read.table(file="Tables/SpParams.txt", sep="\t", header = TRUE)
results_sureau <- read.table(file="Tables/results_sureau.txt", sep="\t", header = TRUE)
results_sperry <- read.table(file="Tables/results_sperry.txt", sep="\t", header = TRUE)
results_sperry_segmented <- read.table(file="Tables/results_sperry_segmented.txt", sep="\t", header = TRUE)

results_sureau <- results_sureau[results_sureau$day_failure_leaf>=50, ]
results_sperry <- results_sperry[results_sperry$day_failure_leaf>=50, ]

results_sureau$model <- "Sureau"
results_sperry$model <- "Sperry"
results <- dplyr::bind_rows(results_sperry, results_sureau)
results <- results |>
  dplyr::rename(Name = Species) |>
  dplyr::left_join(SpParams, by= "Name")

ggplot(SpParams)+
  geom_point(aes(x = VCstem_P50, y = Gs_P50))+
  geom_abline(intercept= 0, slope = 1)+
  theme_bw()


ggplot(results)+
  geom_point(aes(x = VCstem_P50, y = day_closure, col = model))+
  geom_smooth(aes(x = VCstem_P50, y = day_closure, col = model))+
  ylab("Stem P50")+ xlab("Day of stomatal closure")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = VCstem_P50, y = day_failure_stem, col = model))+
  geom_smooth(aes(x = VCstem_P50, y = day_failure_stem, col = model))+
  xlab("Stem P50")+ ylab("Day of hydraulic failure")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.85,0.85))


ggplot(results)+
  geom_point(aes(y = E_max, x = actual_kplant, col = model))+
  xlab("Whole-plant conductance")+ ylab("Max. E (mmol H2O·m-2·s-1)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.85,0.2))

ggplot(results)+
  geom_point(aes(x = E_max, y = day_closure, col = model))+
  ylab("Day of stomatal closure")+ xlab("Max. E (mmol H2O·m-2·s-1)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = actual_kplant, y = day_closure, col = model))+
  ylab("Day of stomatal closure")+ xlab("Whole-plant conductance")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = E_max, y = day_failure_stem, col = model))+
  ylab("Day of hydraulic failure")+ xlab("Max. E (mmol H2O·m-2·s-1)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = actual_kplant, y = day_failure_stem, col = model))+
  ylab("Day of hydraulic failure")+ xlab("Whole-plant conductance")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = psi_min, y = day_failure_stem, col = model))+
  ylab("Day of hydraulic failure")+ xlab("Midday water potential (initial)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = psi_min, y = day_closure, col = model))+
  ylab("Day of stomatal closure")+ xlab("Midday water potential (initial)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.8))

ggplot(results)+
  geom_point(aes(x = VCstem_P50, y = psi_min, col = model))+
  geom_abline(intercept = 0, slope = 1)+
  xlab("Water potential corresponding to 50% stem PLC (MPa)")+ ylab("Initial midday water potential (MPa)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.2))

ggplot(results)+
  geom_point(aes(x = Gs_P50, y = psi_min, col = model))+
  geom_abline(intercept = 0, slope = 1)+
  xlab("Water potential corresponding to 50% Gs (MPa)")+ ylab("Initial midday water potential (MPa)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.2))

ggplot(results)+
  geom_point(aes(x = Nleaf, y = psi_min, col = model))+
  xlab("Leaf nitrogen per dry mass")+ ylab("Initial midday water potential (MPa)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.2))

ggplot(results)+
  geom_point(aes(x = E_max, y = psi_min, col = model))+
  ylab("Midday water potential (initial)")+ xlab("Max. E (mmol H2O·m-2·s-1)")+
  scale_color_discrete("")+
  theme_bw()+theme(legend.position = c(0.8,0.7))





# ggplot()+
#   geom_point(aes(x = results_sureau$day_failure_leaf, y = results_sperry_segmented$day_failure_leaf))+
#   geom_abline(intercept= 0, slope = 1)+
#   ylim(c(0,365))+xlim(c(0,365))+
#   theme_bw()

ggplot()+
  geom_point(aes(x = results_sureau$day_failure_leaf, y = results_sperry$day_failure_leaf))+
  geom_abline(intercept= 0, slope = 1)+
  ylim(c(0,365))+xlim(c(0,365))+
  theme_bw()

results$kplantConsensus <- NewParams$kplantConsensus
results$GS_SM <- NewParams$Gs_P50 - NewParams$P50_VC 
results$P50_VC <- NewParams$P50_VC
results$Gs_P50 <- NewParams$Gs_P50
results$Gswmax <- NewParams$gs_consensus
results$GrowthForm <- NewParams$GrowthForm

results_def <- results[results$Gswmax < 1,]


ggplot(results_def)+
  geom_point(aes(x = kplantConsensus, y = Gswmax))+
  theme_bw()
ggplot(results_def)+
  geom_point(aes(x = actual_kplant_sureau, y = Gswmax))+
  theme_bw()

ggplot(results_def)+
  geom_point(aes(x = P50_VC, y = Gswmax))+
  theme_bw()


ggplot(results_sperry)+
  geom_point(aes(x = actual_kplant, y = Gswmax))+
  theme_bw()


ggplot(results_def)+
  geom_point(aes(x = psi_min_sureau, y = psi_min_sperry, col = P50_VC))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

ggplot(results_def)+
  geom_point(aes(x = E_max_sureau, y = E_max_sperry, col = P50_VC))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

ggplot(results_def)+
  geom_point(aes(x = day_closure_sureau, y = day_closure_sperry, col = P50_VC))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

g4 <- ggplot(results_def)+
  geom_point(aes(x = day_failure_leaf_sureau, y = day_failure_leaf_sperry, col = P50_VC))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()


ggplot(results)+
  geom_point(aes(x = P50_VC, y = day_failure_leaf_sureau))+
  theme_bw()

ggplot(results)+
  geom_point(aes(x = P50_VC, y = day_failure_leaf_sperry))+
  theme_bw()
ggplot(results)+
  geom_point(aes(x = P50_VC, y = day_failure_stem_sperry))+
  theme_bw()

ggplot(results_def)+
  geom_point(aes(x = GS_SM, y = day_failure_leaf_sureau, col = GrowthForm))+
  theme_bw()