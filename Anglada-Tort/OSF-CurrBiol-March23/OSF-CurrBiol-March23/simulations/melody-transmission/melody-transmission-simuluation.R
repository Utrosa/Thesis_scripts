################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis simulations iterated tranmission process
################################################################################


################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(egg)
library(ggpubr)
library(cowplot)
library(viridis)
library(philentropy)  # JSD

# methods
source("utils/plots.R")  # methods for plotting
source("utils/stats.R")  # general methods
source("utils/features.R")  # methods for calculating melodic features

# global parameters
set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())


MAX_INTERVAL_SIZE = 40
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

current_path = "simulations/melody-transmission/"


################################################################################
# Model: interval size
################################################################################
NBOOT = 100
BW = 0.25

data <- read_csv("simulations/melody-transmission/melody-transmisison-1SD-rejection.csv")  %>% 
  separate(target_pitches, c("p1", "p2", "p3"), sep = ",") %>% 
  mutate(p1 = parse_number(p1), p2 = parse_number(p2), p3 = parse_number(p3)) %>% 
  mutate(int1 = p2-p1, int2 = p3-p2) %>% 
  mutate(degree=iteration)
  
 
dat = data %>% 
  # filter(degree %in% 0:9) %>% 
  mutate(
    gg = case_when(
      degree == 0 ~ "Generation 0",
      degree < 10 ~ "10 generations",
      degree < 20 ~ "Generation 11-20",
      degree < 30 ~ "30 generations",
      degree < 40 ~ "Generation 31-40",
      degree < 50 ~ "50 generations",
      degree < 60 ~ "Generation 51-60",
      degree < 70 ~ "70 generations",
      degree < 80 ~ "Generation 71-80",
      degree < 90 ~ "90 generations",
      degree < 100 ~ "Generation 91-100"
    )
  ) %>% 
  filter(gg != "10 generations") %>% 
  filter(gg != "30 generations") %>% 
  filter(gg != "50 generations") %>% 
  filter(gg != "70 generations") %>% 
  filter(gg != "90 generations") %>% 
  pivot_longer(int1:int2)  

# reorder levels
dat$gg = factor(dat$gg, levels = c(
  "Generation 0",
  "Generation 11-20",
  "Generation 31-40",
  "Generation 51-60",
  "Generation 71-80",
  "Generation 91-100"
))

marginals = dat %>%  
  ggplot(aes_string(x="value")) + 
  geom_density(bw=BW)  +
  scale_x_continuous(breaks=seq(min(vertical.lines), max(vertical.lines), 5),
                     limits = c(min(vertical.lines), max(vertical.lines)))  + 
  geom_vline(xintercept = seq(-40,40,5), colour = "lightgrey", linetype="dashed") + 
  facet_wrap(~ gg, ncol = 1) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 12))  

ggsave(paste0(current_path, "simulations_marginals_trial.pdf"),  
       marginals, 
       height = 25, width = 20, units = "cm", dpi=300)


data_ready = data %>% 
  mutate(raw_diff1 = int1 - lag(int1), raw_diff2 = int2 - lag(int2))  %>% 
  unite(raw_diffs, raw_diff1, raw_diff2, sep = ",") %>% 
  mutate(raw_diffs =  map(strsplit(raw_diffs, ","), ~ as.numeric(.))) %>% 
  mutate(rmse_interval = map_dbl(raw_diffs, compute_error))


# run features (takes some time)
features_data_ready = get_bootstrapped_features_simulations(
  data_ready, 
  c("int1", "int2"), 
  NBOOT
)


error = summarise_feature(features_data_ready, RMSE_interval)
int.size = summarise_feature(features_data_ready, abs_int_size)
entropy = summarise_feature(features_data_ready, entropy_0.25)


plot_error = make_lineplot_ribon_features(
  error, degree, mean, sd, "RMSE interval") 

plot_interval = make_lineplot_ribon_features(
  int.size, degree, mean, sd, "Mean Abs Interval Size") 

plot_entropy = make_lineplot_ribon_features(
  entropy, degree, mean, sd, "Interval Entropy") 

plots = plot_grid(plot_entropy, plot_interval, plot_error, ncol = 3)

ggsave(paste0(current_path, "simulations_trial.pdf"),  
       plots, 
       height = 10, width = 25, units = "cm", dpi=300)




