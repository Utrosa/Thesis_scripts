################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis simulations
################################################################################

# Models simulating production and perception (Experiment 4 and 5)

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
library(readxl)

# methods
source("utils/plots.R")  # methods for plotting
source("utils/stats.R")  # general methods

set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

MAX_INTERVAL_SIZE = 13
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

current_path = "simulations/production-perception/"


################################################################################
# Model: interval size
################################################################################
# models
boot_1000_M_intrval <- read_csv("simulations/production-perception/models/simulations-runs.boot.1000.btw-2tones-Interval-size-model.csv")
boot_1000_sing<- read_csv("simulations/production-perception/models/simulations-runs.boot.1000.btw-2tones.csv")
# combined
boot_1000_M_combined <- read_csv("simulations/production-perception/combined/simulations-runs.boot.1000.btw-2tones-Combined-model.csv")


colors_slider = c("grey", "#E13342FF")

sum_1000_M_intrval = prepare_data_sim(boot_1000_M_intrval)
sum_1000_sing = prepare_data_sim(boot_1000_sing)

plot_motor = plot_simulation_models(sum_1000_M_intrval, sum_1000_sing, colors_slider)


ggsave(paste0(current_path, "plot_motor_simulation.pdf"), 
       plot_motor, 
       width = 15, height = 8, units = "cm", dpi=300, bg = "white")


# over generations (SI)
degree_seq = seq(0,10,1)

boot_1000_M_interval_evolution <- read_csv("simulations/production-perception/models/simulations-sample.1000.btw-2tones-Interval-size-model.csv") %>% 
  mutate(degree = degree_seq)
sum_1000_sing_evolution <- read_csv("simulations/production-perception/models/simulations-sample.1000.btw-2tones.csv") %>% 
  mutate(degree = degree_seq)

sum_1000_M_interval_evolution = prepare_data_sim_generations(boot_1000_M_interval_evolution)
sum_1000_sing_evolution = prepare_data_sim_generations(sum_1000_sing_evolution)

plot_motor_evolution = plot_simulation_models_generations(sum_1000_M_interval_evolution, sum_1000_sing_evolution, colors_slider)

ggsave(paste0(current_path, "plot_motor_simulation_generations.pdf"), 
       plot_motor_evolution, 
       width = 14, height = 20, units = "cm", dpi=300, bg = "white")


# fitting function
fun_interval_model <- read_excel("simulations/production-perception/models/poly_fit_interval-size-model.xlsx") %>% 
  summarise_all(mean) %>% 
  pivot_longer(`-14.9`:`14.9`, names_to = "x", values_to = "mean") %>% 
  mutate(x = as.numeric(x))

plot_function_interval = plot_simulation_function(fun_interval_model, colors_slider)

ggsave(paste0(current_path, "plot_function_interval.pdf"), 
       plot_function_interval, 
       width = 14, height = 8, units = "cm", dpi=300, bg = "white")


################################################################################
# Model: preferences
################################################################################
boot_1000_P_preference <- read_csv("simulations/production-perception/models/simulations-runs.boot.1000.btw-2tones-Preference-model.csv")

colors_slider = c("grey", "#35B779FF")

sum_1000_P_preference = prepare_data_sim(boot_1000_P_preference)

plot_preference = plot_simulation_models(sum_1000_P_preference, sum_1000_sing, colors_slider)


ggsave(paste0(current_path, "plot_preference_simulation.pdf"), 
       plot_preference, 
       width = 15, height = 8, units = "cm", dpi=300, bg = "white")


# over generations (SI)
boot_1000_M_preference <- read_csv("simulations/production-perception/models/simulations-sample.1000.btw-2tones-Preference-model.csv") %>% 
  mutate(degree = degree_seq)
sum_1000_sing_evolution <- read_csv("simulations/production-perception/models/simulations-sample.1000.btw-2tones.csv") %>% 
  mutate(degree = degree_seq)

sum_1000_M_preferemce_evolution = prepare_data_sim_generations(boot_1000_M_preference)
sum_1000_sing_evolution = prepare_data_sim_generations(sum_1000_sing_evolution)

plot_preference_evolution = plot_simulation_models_generations(sum_1000_M_preferemce_evolution, sum_1000_sing_evolution, colors_slider)

ggsave(paste0(current_path, "plot_preference_simulation_generations.pdf"), 
       plot_preference_evolution, 
       width = 14, height = 20, units = "cm", dpi=300, bg = "white")


# fitting function
fun_preference_model <- read_excel("simulations/production-perception/models/pref_smoothes_preference-model.xlsx") %>% 
  summarise_all(mean) %>% 
  pivot_longer(`-14.9`:`14.9`, names_to = "x", values_to = "mean") %>% 
  mutate(x = as.numeric(x))

plot_function_preference = plot_simulation_function(fun_preference_model, colors_slider)

ggsave(paste0(current_path, "plot_function_preference.pdf"), 
       plot_function_preference, 
       width = 14, height = 8, units = "cm", dpi=300, bg = "white")


################################################################################
# Combined
################################################################################
boot_1000_Combined <- read_csv("simulations/production-perception/combined/simulations-runs.boot.1000.btw-2tones-Combined-model.csv")


# colors
library(scales)
show_col(viridis_pal(option = "C")(10))
colors_combined = c("grey", "#FA9E3BFF")

sum_1000_Combined = prepare_data_sim(boot_1000_Combined)

plot_combined = plot_simulation_models(sum_1000_Combined, sum_1000_sing, colors_combined)

ggsave(paste0(current_path, "plot_combined_simulation.pdf"), 
       plot_combined, 
       width = 15, height = 8, units = "cm", dpi=300, bg = "white")


# over generations (SI)
boot_1000_Combined <- read_csv("simulations/production-perception/combined/simulations-sample.1000.btw-2tones-Combined-model.csv") %>% 
  mutate(degree = degree_seq)
sum_1000_sing_evolution <- read_csv("simulations/production-perception/combined/simulations-sample.1000.btw-2tones.csv") %>% 
  mutate(degree = degree_seq)

sum_1000_M_combined_evolution = prepare_data_sim_generations(boot_1000_Combined)
sum_1000_sing_evolution = prepare_data_sim_generations(sum_1000_sing_evolution)

plot_preference_evolution = plot_simulation_models_generations(
  sum_1000_M_combined_evolution, sum_1000_sing_evolution, colors_combined)

ggsave(paste0(current_path, "plot_combined_simulation_generations.pdf"), 
       plot_preference_evolution, 
       width = 14, height = 20, units = "cm", dpi=300, bg = "white")


