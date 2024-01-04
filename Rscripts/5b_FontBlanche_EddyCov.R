# Update medfate first from devel branch!!!
# remotes::install_github("emf-creaf/medfate", ref = "devel")
library(medfate)
library(readr)
library(tidyverse)
library(cowplot)

source("Rscripts/0b_Ancillary_FontBlanche.R")


# Terrain -----------------------------------------------------------------
fb_latitude <- 43.24
fb_elevation <- 420
fb_slope <- 0 
fb_aspect <- 0

# Meteo -------------------------------------------------------------------
fb_meteo_raw <- read_delim("Data/FontBlanche/Climate_FontBlanche_GapFilled.csv",
                       delim = ";", escape_double = FALSE, trim_ws = TRUE)
fb_meteo <- read_delim("Data/FontBlanche/Climate_FontBlanche_GapFilled.csv",
                       delim = ";", escape_double = FALSE, trim_ws = TRUE) |>
  rename(dates = Date,
         MinTemperature = Tair_min,
         MaxTemperature = Tair_max,
         MeanTemperature = Tair_mean,
         Radiation = RG_sum,
         Precipitation = PPT_sum, 
         MinRelativeHumidity = RHair_min, 
         MaxRelativeHumidity = RHair_max, 
         MeanRelativeHumidity = RHair_mean,
         WindSpeed = WS_mean) |>
  mutate(dates = as.Date(dates, format = "%d/%m/%Y"),
         Radiation = Radiation/100)|>
  filter(dates >= as.Date("2016-01-01"),
         dates <= as.Date("2022-12-31"))


# Cochard initialization and run ----------------------------------------------------------
#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefillStem <- "annual"
control$cavitationRefillLeaves <- "annual"
control$bareSoilEvaporation <- TRUE
control$plantCapacitance <- TRUE
control$cavitationFlux <- TRUE
control$sapFluidityVariation <- TRUE
control$leafCuticularTranspiration <- TRUE
control$stemCuticularTranspiration <- TRUE
control$sunlitShade <- TRUE
control$gs_NightFrac <- 0.001
control$stomatalSubmodel <- "Baldocchi"
control$rhizosphereOverlap <- "partial"
control$ndailysteps <- 48 # Half hour time steps
x2p <- fontblanche_input(control)
S2p <- spwb(x2p, fb_meteo, 
            latitude = fb_latitude, elevation = fb_elevation,
            slope = fb_slope, aspect = fb_aspect)
saveRDS(S2p, "Rdata/FontBlanche/FontBlanche_Sureau_EddyCov.rds")

# Eddy Covariance fluxes --------------------------------------------------
EC_FB <- read_csv("Data/FontBlanche/FR-FBn_2016_2022.csv", 
                             col_types = cols("...1" = col_skip(),
                                              "station.name" = col_skip()))


# Extract subdaily canopy energy balance ----------------------------------
CEB <- medfate::extractSubdaily(S2p, "CanopyEnergyBalance")|>
  tibble::as_tibble()

LE_30min <- data.frame(datetime = EC_FB$DateTime, EC = EC_FB$LE_F, PRED = CEB$LEcan)
LE_day <- LE_30min |>
  mutate(date = as.Date(datetime)) |>
  group_by(date) |>
  summarise(EC = mean(EC), PRED = mean(PRED)) |>
  mutate(RainyDay = fb_meteo$Precipitation>0,
         EC_NR = EC,
         PRED_NR = PRED)

LE_day$EC_NR[LE_day$RainyDay] <- NA
LE_day$PRED_NR[LE_day$RainyDay] <- NA

g_scatter <- ggplot(LE_day, aes(x = EC, y = PRED))+
  geom_point(aes(col= RainyDay), size = 0.5)+
  geom_abline(col = "black")+
  xlim(c(0,120))+ ylim(c(0,120))+
  xlab("Observed LE [W/s]")+ylab("Predicted LE [W/s]")+
  scale_color_manual("", labels = c("Dry day", "Wet day"), values=c("red", "blue"))+
  theme_bw()+ theme(legend.position = c(0.85,0.10))
ggsave2("Plots/FontBlanche_LE_scatter.png", g_scatter, width = 5, height=5)

g_dyn <- LE_day |>
  select(date, EC, PRED)|>
  pivot_longer(cols = c("EC", "PRED"), 
               names_to = "method", values_to = "LE") |>
  ggplot()+
  geom_line(aes(x = date, y = LE, col=method), alpha = 0.5)+
  xlab("")+ylab("LE [W/s]")+
  scale_color_manual(NULL, labels = c("Observed", "Predicted"), 
                     values = c("black", "red"))+
  theme_bw()+ theme(legend.position = c(0.90,0.87))
ggsave2("Plots/FontBlanche_LE_dynamics.png", g_dyn, width = 8, height=4)

g_dyn_NR <- LE_day |>
  select(date, EC_NR, PRED_NR)|>
  pivot_longer(cols = c("EC_NR", "PRED_NR"), 
               names_to = "method", values_to = "LE") |>
  ggplot()+
  geom_line(aes(x = date, y = LE, col=method), alpha = 0.5)+
  xlab("")+ylab("LE [W/s]")+
  scale_color_manual(NULL, labels = c("Observed", "Predicted"), 
                     values = c("black", "red"))+
  theme_bw()+ theme(legend.position = c(0.90,0.87))
ggsave2("Plots/FontBlanche_LE_dynamics_NoRainy.png", g_dyn_NR, width = 8, height=4)

stats<- function(observed, predicted) {
  is_mis <- is.na(observed) | is.na(predicted)
  observed <- observed[!is_mis]
  predicted <- predicted[!is_mis]
  list(
    bias = mean(predicted - observed, na.rm=TRUE),
    rmse = sqrt(mean((predicted - observed)^2, na.rm=TRUE)),
    r2 = cor(predicted, observed, use = "complete")^2
  )
}

stats(LE_day$EC, LE_day$PRED)
stats(LE_day$EC_NR, LE_day$PRED_NR)
