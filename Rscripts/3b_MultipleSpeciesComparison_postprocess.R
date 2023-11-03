library(medfate)
library(ggplot2)

# Process results ---------------------------------------------------------
SpParams <- read.table(file="Tables/SpParams.txt", sep="\t", header = TRUE)
results_sureau <- read.table(file="Tables/results_sureau.txt", sep="\t", header = TRUE)
results_sperry <- read.table(file="Tables/results_sperry.txt", sep="\t", header = TRUE)
results_sperry_segmented <- read.table(file="Tables/results_sperry_segmented.txt", sep="\t", header = TRUE)



ggplot(SpParams)+
  geom_point(aes(x = VCstem_P50, y = Gs_P50))+
  geom_abline(intercept= 0, slope = 1)+
  theme_bw()

results_sperry <- results_sperry |>
  dplyr::rename(Name = Species) |>
  dplyr::left_join(SpParams, by= "Name")

results_sperry_segmented <- results_sperry_segmented |>
  dplyr::rename(Name = Species) |>
  dplyr::left_join(SpParams, by= "Name")

results_sureau <- results_sureau |>
  dplyr::rename(Name = Species) |>
  dplyr::left_join(SpParams, by= "Name")

ggplot(results_sperry)+
  geom_point(aes(x = actual_kplant, y = E_max))+
  theme_bw()

ggplot(results_sperry)+
  geom_point(aes(x = E_max, y = day_closure))+
  theme_bw()

ggplot(results_sureau)+
  geom_point(aes(x = E_max, y = day_closure))+
  theme_bw()

ggplot(results_sperry_segmented)+
  geom_point(aes(x = E_max, y = day_closure))+
  theme_bw()


ggplot()+
  geom_point(aes(x = results_sureau$day_closure, y = results_sperry$day_closure), col = "black")+
  geom_point(aes(x = results_sureau$day_closure, y = results_sperry_segmented$day_closure), col = "red")+
  geom_abline(intercept= 0, slope = 1)+
  ylim(c(0,365))+xlim(c(0,365))+
  theme_bw()

ggplot()+
  geom_point(aes(x = results_sureau$day_failure_leaf, y = results_sperry$day_failure_leaf), col = "black")+
  geom_point(aes(x = results_sureau$day_failure_leaf, y = results_sperry_segmented$day_failure_leaf), col = "red")+
  geom_abline(intercept= 0, slope = 1)+
  ylim(c(0,365))+xlim(c(0,365))+
  theme_bw()

ggplot()+
  geom_point(aes(x = results_sureau$day_failure_leaf, y = results_sperry_segmented$day_failure_leaf))+
  geom_abline(intercept= 0, slope = 1)+
  ylim(c(0,365))+xlim(c(0,365))+
  theme_bw()

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