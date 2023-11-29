library(tidyverse)
library(ggplot2)
library(cowplot)

model = "SurEau"

sa_timetofailure <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_timetofailure.rds"))
sa_timetoclosure <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_timetoclosure.rds"))
sa_survivaltime <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_survivaltime.rds"))
sa_gpp <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_gpp.rds"))

sa_timetofailure_red <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_timetofailure_red.rds"))
sa_timetoclosure_red <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_timetoclosure_red.rds"))
sa_survivaltime_red <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_survivaltime_red.rds"))
sa_gpp_red <- readRDS(paste0("Rdata/sensitivity/", model, "/sa_salt_gpp_red.rds"))

if(model=="SurEau") {
  parNames = c("TAW",
               "LAI", 
               "Z95", 
               "Plant_kmax", 
               "Gswmax", 
               "Vmax298/Jmax298", 
               "Gs_P50", 
               "Gswmin", 
               "VC_P50",
               "Vsapwood")
  red_exclude <- c(1,2)
  parNames_red <- parNames[-red_exclude]
} else {
  parNames = c("TAW",
               "LAI", 
               "Z95", 
               "Plant_kmax", 
               "Gswmax", 
               "Vmax298/Jmax298", 
               "Gswmin", 
               "VC_d")
  red_exclude <- c(1,2)
  parNames_red <- parNames[-red_exclude]
}

sens_tidy <- function(x, response) {
  data.frame(ParName = parNames, Response = response, Total = x$T$original, FirstOrder = x$S$original) 
}
sens_tidy_red <- function(x, response) {
  data.frame(ParName = parNames_red, Response = response, Total = x$T$original, FirstOrder = x$S$original) 
}

ttc_df <- sens_tidy(sa_timetoclosure, response = "TSC")
thf_df <- sens_tidy(sa_timetofailure, response = "THF")
st_df <- sens_tidy(sa_survivaltime, response = "ST")
gpp_df <- sens_tidy(sa_gpp, response = "GPP")
all_df <- bind_rows(ttc_df, thf_df, st_df, gpp_df) |>
  mutate(Response = factor(Response, levels = c("GPP", "TSC", "THF", "ST")),
         ParName = factor(ParName, levels = parNames))

ttc_red_df <- sens_tidy_red(sa_timetoclosure_red, response = "TSC")
thf_red_df <- sens_tidy_red(sa_timetofailure_red, response = "THF")
st_red_df <- sens_tidy_red(sa_survivaltime_red, response = "ST")
gpp_red_df <- sens_tidy_red(sa_gpp_red, response = "GPP")
all_red_df <- bind_rows(ttc_red_df, thf_red_df, st_red_df, gpp_red_df) |>
  mutate(Response = factor(Response, levels = c("GPP","TSC", "THF", "ST")),
         ParName = factor(ParName, levels = parNames_red))


p1<-ggplot(all_df) +
  geom_bar(aes(x = ParName, y = Total, fill = Response), stat = "identity", position="dodge")+
  coord_flip()+
  scale_x_discrete(limits = parNames[length(parNames):1])+
  scale_fill_brewer(type = "qual")+
  ylab("Total effect")+ xlab("")+labs(title = "Total effect")+
  theme_bw()

p2<-ggplot(all_df) +
  geom_bar(aes(x = ParName, y = FirstOrder, fill = Response), stat = "identity", position="dodge")+
  coord_flip()+
  scale_x_discrete(limits = parNames[length(parNames):1])+
  scale_fill_brewer(type = "qual")+
  ylab("First-order effect")+ xlab("")+labs(title = "First-order effect")+
  theme_bw()
p<-plot_grid(p1+theme(legend.position = "none"), 
          p2+theme(legend.position = "none"), 
          get_legend(p1), nrow = 1, rel_widths = c(1,1,0.2))
ggsave2(paste0("Plots/Sensitivity_", model, "_all.png"), p, width = 12, height= 7, bg = "white")

p3<- ggplot(all_red_df) +
  geom_bar(aes(x = ParName, y = Total, fill = Response), stat ="identity", position="dodge")+
  coord_flip()+
  scale_x_discrete(limits = parNames_red[length(parNames_red):1])+
  scale_fill_brewer(type = "qual")+
  ylab("Total effect")+xlab("")+labs(title = "Total effect")+
  theme_bw()
p4<-ggplot(all_red_df) +
  geom_bar(aes(x = ParName, y = FirstOrder, fill = Response), stat ="identity", position="dodge")+
  coord_flip()+
  scale_x_discrete(limits = parNames_red[length(parNames_red):1])+
  scale_fill_brewer(type = "qual")+
  ylab("First-order effect")+xlab("")+labs(title = "First-order effect")+
  theme_bw()
p<-plot_grid(p3+theme(legend.position = "none"), 
          p4+theme(legend.position = "none"), 
          get_legend(p3), nrow = 1, rel_widths = c(1,1,0.2))
ggsave2(paste0("Plots/Sensitivity_", model, "_red.png"), p, width = 12, height= 6, bg = "white")
