################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis cross-cultural comparison
################################################################################

# Structure: (1) distributions, (2) features, (3) stats
# Experiments 1 and 10

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

data_exp10 <- read_csv("data/experiment10/exp10_btw_3tones_india_120chains.csv") %>% 
  group_by(network_id) %>% 
  mutate(n_trials = n()) %>% 
  filter(n_trials > 9)  # completing all chains in this experiment was not possible because the limited pool of active participants

length(unique(data_exp10$participant_id)) # 54
length(unique(data_exp10$network_id)) # 120
table(data_exp10$degree) # 95

current_path = "analysis/5-cross-cultural/"

# colors
library(scales)
show_col(viridis_pal(option = "D")(16))
#  groups
colors_culture = c("grey", "#54C568FF")



# peaks from boot peakfinding method
peaks_exp1_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp10_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-india-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

peaks_exp1 = get_sig_peaks(peaks_exp1_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks_exp10 = get_sig_peaks(peaks_exp10_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))


# individual peaks (interval 1 and 2):
# peaks from boot peakfinding method
peaks_exp10_int1 = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-india-interval-1.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp10_int2 = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-india-interval-2.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))


################################################################################
# Distribution
################################################################################
MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)
sung_intervals = c("sung_interval1", "sung_interval2")
BW = 0.25


# 2D KDE
peaks1 = get_sig_peaks(peaks_exp10_int1, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks2 = get_sig_peaks(peaks_exp10_int2, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))

kde_plot_main <- generate_2Dkde_main(data_exp10, peaks1, peaks2, interval_range, BW)

ggsave(paste0(current_path, "2DKDE_main.pdf"), 
       kde_plot_main, 
       height = 18, width = 18,
       units = "cm",
       dpi=300)


# 2KDEs
plot_2kde_us_basic = data_exp1 %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")

plot_2kde_india_basic = data_exp10 %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")

ggsave(paste0(current_path, "plot_2kde_us_basic.pdf"), 
       plot_2kde_us_basic, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)

ggsave(paste0(current_path, "plot_2kde_india_basic.pdf"), 
       plot_2kde_india_basic, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)

# Comparing 1KDE separately
joint.marginals_us = make_marginals_1datasets(
  data_exp1, 
  data_exp10,
  peaks_exp1,
  sung_intervals,
  NBOOT, BW, interval_range,
  colors_culture
  ) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.05), limits = c(0, 0.17))

joint.marginals_india = make_marginals_1datasets(
  data_exp10, 
  data_exp1,
  peaks_exp10,
  sung_intervals,
  NBOOT, BW, interval_range,
  c("#54C568FF", "black")
) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.05), limits = c(0, 0.17))


ggsave(paste0(current_path, "joint.marginals_us.pdf"), 
       joint.marginals_us, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")
ggsave(paste0(current_path, "joint.marginals_india.pdf"), 
       joint.marginals_india, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")


################################################################################
# Features
################################################################################
# peaks from matlab procedure
path_peaks = "peaks/number_peaks_iterations-btw-3tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp1 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-btw-3tones-india-all-intervals.NBOOT-1000.txt"
boot_peaks_exp10 = read_boot_peaks_data(path_peaks) 

# features
features_exp1 <- read_csv("features/exp1_features_btw_3tones_590chians_1kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp1)

features_exp10 <- read_csv("features/exp10_features_btw_3tones_india_120chains_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp10)


features_exp1$app = "aUS"
features_exp10$app = "India"

features_compare = as_tibble(rbind(features_exp10, features_exp1))

vars_to_sum = c(
  "RMSE_interval",
  "abs_int_size",
  "n_peaks",
  "entropy_0.25"
)


plots_features_compare = plot_features_2datasets(
  features_compare, 
  vars_to_sum,
  normalize = FALSE,
  colors_culture, 
  is_viridis = FALSE,
  legend_off = TRUE
)

plots_features_compare$plot_entropy_0.25
plots_features_compare$plot_peaks
plots_features_compare$plot_int.size

plots_features_compare$plot_error

# normalize (normalize = TRUE)
ggsave(paste0(current_path, "final_entropy.pdf"), 
       plots_features_compare$plot_entropy_0.25, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "final_peaks.pdf"), 
       plots_features_compare$plot_peaks,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "final_interval_size.pdf"), 
       plots_features_compare$plot_int.size,
       height = 10, width = 10, units = "cm", dpi=300)

# do not normalize (normalize = FALSE)
ggsave(paste0(current_path, "final_error.pdf"), 
       plots_features_compare$plot_error,
       height = 10, width = 10, units = "cm", dpi=300)


# mean interval size in last generation
features_compare %>% 
  filter(degree == 10) %>% 
  group_by(app) %>% 
  summarise(
    m = mean(abs_int_size),
    sd = sd(abs_int_size) * 1.96,
    low_ci = m - sd,
    high_ci = m + sd)


################################################################################
# stats
################################################################################
boot_stats_exp10 = run_lm_bootraped_data(features_exp10)

# stats 
boot_stats_exp10[[2]]

