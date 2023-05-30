################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis short-melodies (Experiment 1)
################################################################################

# Structure: (1) distributions, (2) features, (3) stats
# Experiment 1

################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(egg)
library(ggpubr)
library(cowplot)
library(viridis)

# methods
source("utils/plots.R")  # methods for plotting
source("utils/stats.R")  # general methods

# global parameters
set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

# data
data_exp1 <- read_csv("data/experiment1/exp1_btw_3tones_590chains.csv")

# peaks from boot peakfinding method
peaks_exp1_int1 = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-interval-1.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp1_int2 = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-interval-2.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

current_path = "analysis/1-short-melodies/"

################################################################################
# 1. Distributions
################################################################################
# 2D KDE
BW = 0.25

peaks1 = get_sig_peaks(peaks_exp1_int1, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks2 = get_sig_peaks(peaks_exp1_int2, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))

kde_plot_main <- generate_2Dkde_main(data_exp1, peaks1, peaks2, interval_range, BW)

ggsave(paste0(current_path, "2DKDE_main.pdf"), 
       kde_plot_main, 
       height = 18, width = 18,
       units = "cm",
       dpi=300)


# marginals over generations
kde_plot_evolution_btw = generate_2Dkde_facet(data_exp1, interval_range)

ggsave(paste0(current_path, "2DKDE_evolution_btw_revised.pdf"), 
       kde_plot_evolution_btw, 
       height = 15, width = 30,
       units = "cm",
       dpi=300)


################################################################################
# 2. Features
################################################################################
# peaks from matlab procedure
path_peaks = "peaks/number_peaks_iterations-btw-3tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp1 = read_boot_peaks_data(path_peaks) 

# melodic features
features_exp1 <- read_csv("features/exp1_features_btw_3tones_590chians_1kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp1)


# plot features
plot_features_exp1 = plot_features(features_exp1)

plot_features_exp1$plot_entropy_0.25
plot_features_exp1$plot_peaks
plot_features_exp1$plot_int.size
plot_features_exp1$plot_error


ggsave(paste0(current_path, "entropy_0.25.pdf"),  
       plot_features_exp1$plot_entropy_0.25, 
       height = 8, width = 8, units = "cm", dpi=300)
ggsave(paste0(current_path, "peaks.pdf"),  
       plot_features_exp1$plot_peaks,
       height = 8, width = 8, units = "cm", dpi=300)
ggsave(paste0(current_path, "interval_size.pdf"),  
       plot_features_exp1$plot_int.size,
       height = 8, width = 8, units = "cm", dpi=300)
ggsave(paste0(current_path, "error.pdf"), 
       plot_features_exp1$plot_error, 
       height = 8, width = 8, units = "cm", dpi=300)


################################################################################
# 3. Stats
################################################################################

boot_stats_exp1 = run_lm_bootraped_data(features_exp1)

# stats 
boot_stats_exp1[[2]]




