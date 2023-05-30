################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis memory experiments
################################################################################

# Structure: (1) distributions, (2) features, (3) stats
# Experiments 7-9

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
sung_intervals = c("sung_interval1", "sung_interval2")

# data
data_exp7 <- read_csv("data/experiment7-9/exp7_btw_3tones_memory.5sec_240chains.csv")
data_exp8 <- read_csv("data/experiment7-9/exp8_btw_3tones_memory.10sec_240chains.csv")
data_exp9 <- read_csv("data/experiment7-9/exp9_btw_3tones_memory.control_240chains.csv")

current_path = "analysis/4-memory/"


# peaks from boot peakfinding method
peaks_exp7_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-mem-3tones-5-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp8_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-mem-3tones-10-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

peaks_exp7 = get_sig_peaks(peaks_exp7_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks_exp8 = get_sig_peaks(peaks_exp8_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))


# colors
library(scales)
show_col(viridis_pal(option = "G")(4))

# 3 groups
colors_memory = c("grey", "#38AAACFF", "#40498EFF")

################################################################################
# 1. Distributions
################################################################################

# 2KDE
MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)
BW = 0.25
NBOOT = 1000


# 5 sec
kde_seed_mem5sec = data_exp7 %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "none") 
# 10 sec
kde_seed_mem10sec = data_exp8 %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "none") 
# control
kde_seed_mem.control = data_exp9 %>% 
  filter(degree %in% 7:10) %>% 
  generate_2Dkde_basic(interval_range, "none") 


ggsave(paste0(current_path, "2kde_seed_mem5sec.pdf"), kde_seed_mem5sec, height = 15, width = 15,
       units = "cm", dpi=300)
ggsave(paste0(current_path, "2kde_first3_mem10sec.pdf"), kde_seed_mem10sec, height = 15, width = 15,
       units = "cm", dpi=300)
ggsave(paste0(current_path, "2kde_last3_mem.control.pdf"), kde_seed_mem.control, height = 15, width = 15,
       units = "cm", dpi=300)


################################################################################

# Supplementary

MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

# 5 sec
kde_plot_evolution_5sec = generate_2Dkde_facet_3groups(data_exp7, interval_range) 
# 10 sec
kde_plot_evolution_10sec = generate_2Dkde_facet_3groups(data_exp8, interval_range)
# control
kde_plot_evolution_mem.control = generate_2Dkde_facet_3groups(data_exp9, interval_range)

ggsave(paste0(current_path, "SI_2DKDE_evolution_5sec_revised.pdf"), 
       kde_plot_evolution_5sec, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)

ggsave(paste0(current_path, "SI_2DKDE_evolution_10sec_revised.pdf"), 
       kde_plot_evolution_10sec, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)

ggsave(paste0(current_path, "SI_2DKDE_evolution_mem.control_revised.pdf"), 
       kde_plot_evolution_mem.control, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)


################################################################################
# 2. Features
################################################################################
# peaks from matlab procedure
path_peaks = "peaks/number_peaks_iterations-mem-3tones-5-all-intervals.NBOOT-1000.txt"
boot_peaks_exp7 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-mem-3tones-10-all-intervals.NBOOT-1000.txt"
boot_peaks_exp8 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-mem-3tones-match-all-intervals.NBOOT-1000.txt"
boot_peaks_exp9 = read_boot_peaks_data(path_peaks) 

# features
exp7_features <- read_csv("features/exp7_features_btw_3tones_memory_240chians_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp7)

exp8_features <- read_csv("features/exp8_features_btw_3tones_memory_240chians_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp8)

exp9_features <- read_csv("features/exp9_features_btw_3tones_memory_240chians_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree))%>% 
  # add peaks
  bind_cols(boot_peaks_exp9)

exp7_features$app = "5 sec"
exp8_features$app = "10 sec"
exp9_features$app = "control"

features_compare = as_tibble(rbind(exp7_features, exp8_features, exp9_features))

features_compare$app = factor(features_compare$app,
                                 levels = c("control", "5 sec", "10 sec"))

vars_to_sum = c(
  "RMSE_interval",
  "abs_int_size",
  "n_peaks",
  "entropy_0.25"
)

plots_features_compare = plot_features_2datasets(
  features_compare, 
  vars_to_sum,
  normalize = TRUE,
  colors_memory, 
  is_viridis = FALSE,
  legend_off = TRUE
)

plots_features_compare$plot_entropy_0.25
plots_features_compare$plot_peaks
plots_features_compare$plot_int.size

# normalize (normalize = TRUE)
ggsave(paste0(current_path, "entropy.pdf"), 
       plots_features_compare$plot_entropy_0.25, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "peaks.pdf"), 
       plots_features_compare$plot_peaks,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "interval_size.pdf"), 
       plots_features_compare$plot_int.size,
       height = 10, width = 10, units = "cm", dpi=300)

plots_features_compare$plot_error

# do not normalize (normalize = FALSE)
ggsave(paste0(current_path, "error.pdf"), 
       plots_features_compare$plot_error,
       height = 10, width = 10, units = "cm", dpi=300)


################################################################################
# 3. Stats
################################################################################

boot_stats_exp7 = run_lm_bootraped_data(exp7_features)
boot_stats_exp8 = run_lm_bootraped_data(exp8_features)
boot_stats_exp9 = run_lm_bootraped_data(exp9_features)

# stats 
boot_stats_exp7[[2]]
boot_stats_exp8[[2]]
boot_stats_exp9[[2]]







