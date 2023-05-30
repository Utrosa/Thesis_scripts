################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis social vs asocial transmisison condition in India and the US
################################################################################

# Structure: (1) distributions, (2) features, (3) stats
# Experiment 1 and 10-12

################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(egg)
library(ggpubr)
library(cowplot)
library(viridis)

# Experiments 1 and 10-12

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

# global parameters
set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

# data
data_exp1_us.btw <- read_csv("data/experiment1/exp1_btw_3tones_590chains.csv")
data_exp11_us.wth <- read_csv("data/experiment11/exp11_wth_3tones_615chains.csv")

data_exp10_in.btw <- read_csv("data/experiment10/exp10_btw_3tones_india_120chains.csv") %>% 
  group_by(network_id) %>% 
  mutate(n_trials = n()) %>% 
  filter(n_trials > 9)  # completing all chains in this experiment was not possible because the limited pool of active participants

data_exp12_in.wth <- read_csv("data/experiment12/exp12_wth_3tones_india_223chains.csv") %>% 
  group_by(network_id) %>% 
  mutate(n_trials = n()) %>% 
  filter(n_trials > 10)  # only accept full chains


current_path = "analysis/6-social-individual/"


# peaks from boot peakfinding method
peaks_exp1_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp10_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-btw-3tones-india-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp11_us.wth_raw = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-wth-3tones-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))
peaks_exp12_in.wth_raw = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-wth-3tones-india-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

peaks_exp1_btw = get_sig_peaks(peaks_exp1_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks_exp10_btw = get_sig_peaks(peaks_exp10_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks_exp11_us.wth = get_sig_peaks(peaks_exp11_us.wth_raw, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))
peaks_exp12_in.wth = get_sig_peaks(peaks_exp12_in.wth_raw, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))


################################################################################
# Distribution
################################################################################
MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)
sung_intervals = c("sung_interval1", "sung_interval2")
BW = 0.25

# basic 2KDEs
plot_2kde_exp1_us.btw = data_exp1_us.btw %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")
plot_2kde_exp11_us.wth = data_exp11_us.wth %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")
plot_2kde_exp10_in.btw = data_exp10_in.btw %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")
plot_2kde_exp12_in.wth = data_exp12_in.wth %>% 
  filter(degree %in% 8:10) %>% 
  generate_2Dkde_basic(interval_range, "right")


ggsave(paste0(current_path, "plot_2kde_exp1_us.btw.pdf"), 
       plot_2kde_exp1_us.btw, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "plot_2kde_exp11_us.wth.pdf"), 
       plot_2kde_exp11_us.wth, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "plot_2kde_exp10_in.btw.pdf"), 
       plot_2kde_exp10_in.btw, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "plot_2kde_exp12_in.wth.pdf"), 
       plot_2kde_exp12_in.wth, 
       height = 15, width = 15,
       units = "cm",
       dpi=300)


# marginals
MAX_INTERVAL_SIZE = 13
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)
sung_intervals = c("sung_interval1", "sung_interval2")
BW = 0.25
NBOOT = 1000

# colors
library(scales)
show_col(viridis_pal(option = "B")(8))
#  groups
colors_us = c("#9C179EFF", "grey")
colors_in = c("#E13342FF", "#54C568FF")

colors_wth = c(colors_us[1], colors_in[1])
colors_btw = c(colors_us[2], colors_in[2])


# separated  marginals: WITHIN
joint.marginals_us_wth = make_marginals_1datasets(
  data_exp11_us.wth, 
  data_exp12_in.wth,
  peaks_exp11_us.wth,
  sung_intervals,
  NBOOT, BW, interval_range,
  colors_wth
) +
  scale_y_continuous(breaks=seq(0, 0.1, 0.05), limits = c(0, 0.125))

joint.marginals_india_wth = make_marginals_1datasets(
  data_exp12_in.wth, 
  data_exp11_us.wth,
  peaks_exp12_in.wth,
  sung_intervals,
  NBOOT, BW, interval_range,
  rev(colors_wth)
) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.05), limits = c(0, 0.125))


ggsave(paste0(current_path, "joint.marginals_us_wth.pdf"), 
       joint.marginals_us_wth, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")
ggsave(paste0(current_path, "joint.marginals_india_wth.pdf"), 
       joint.marginals_india_wth, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")


# separated joint marginals: ACROSS
joint.marginals_us = make_marginals_1datasets(
  data_exp1_us.btw, 
  data_exp10_in.btw,
  peaks_exp1_btw,
  sung_intervals,
  NBOOT, BW, interval_range,
  colors_btw
) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.05), limits = c(0, 0.17))

joint.marginals_india = make_marginals_1datasets(
  data_exp10_in.btw, 
  data_exp1_us.btw,
  peaks_exp10_btw,
  sung_intervals,
  NBOOT, BW, interval_range,
  c(colors_btw[2], "black")
) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.05), limits = c(0, 0.17))


ggsave(paste0(current_path, "joint.marginals_us_btw.pdf"), 
       joint.marginals_us, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")
ggsave(paste0(current_path, "joint.marginals_india_btw.pdf"), 
       joint.marginals_india, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")


################################################################################
# supplementary
MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

# marginals over generations vertically 
kde_plot_evolution_btw_3slides_exp1 = generate_2Dkde_facet_3groups(data_exp1_us.btw, interval_range)
kde_plot_evolution_btw_3slides_exp10 = generate_2Dkde_facet_3groups(data_exp10_in.btw, interval_range)
kde_plot_evolution_wth_3slides_exp11 = generate_2Dkde_facet_3groups(data_exp11_us.wth, interval_range)
kde_plot_evolution_wth_3slides_exp12 = generate_2Dkde_facet_3groups(data_exp12_in.wth, interval_range)


ggsave(paste0(current_path, "SI_2DKDE_evolution_3slides_EXP1_US.btw.pdf"), 
       kde_plot_evolution_btw_3slides_exp1, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "SI_2DKDE_evolution_3slides_EXP10_IN.btw.pdf"), 
       kde_plot_evolution_btw_3slides_exp10, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "SI_2DKDE_evolution_3slides_EXP11_US.wth.pdf"), 
       kde_plot_evolution_wth_3slides_exp11, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)
ggsave(paste0(current_path, "SI_2DKDE_evolution_3slides_EXP12_IN.wth.pdf"), 
       kde_plot_evolution_wth_3slides_exp12, 
       height = 30, width = 15,
       units = "cm",
       dpi=300)


################################################################################
# Features
################################################################################
# peaks from matlab procedure
path_peaks = "peaks/number_peaks_iterations-btw-3tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp1 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-btw-3tones-india-all-intervals.NBOOT-1000.txt"
boot_peaks_exp10 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-wth-3tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp11 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-wth-3tones-india-all-intervals.NBOOT-1000.txt"
boot_peaks_exp12 = read_boot_peaks_data(path_peaks) 

# features
features_exp1 <- read_csv("features/exp1_features_btw_3tones_590chians_1kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp1)

features_exp10 <- read_csv("features/exp10_features_btw_3tones_india_120chains_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp10)

features_exp11 <- read_csv("features/exp11_features_wth_3tones_615chians_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp11)

features_exp12 <- read_csv("features/exp12_features_wth_3tones_india_223chains_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp12)

features_exp1$app = "US-across"
features_exp11$app = "US-within"
features_exp10$app = "IN-across"
features_exp12$app = "IN-within"

features_compare_btw = as_tibble(rbind(features_exp1, features_exp10))
features_compare_wth = as_tibble(rbind(features_exp11, features_exp12))

features_compare_us = as_tibble(rbind(features_exp1, features_exp11))
features_compare_india = as_tibble(rbind(features_exp10, features_exp12))

vars_to_sum = c(
  "RMSE_interval",
  "abs_int_size",
  "n_peaks",
  "entropy_0.25"
)


# us
plots_features_compare_us = plot_features_2datasets(
  features_compare_us, 
  vars_to_sum,
  normalize = TRUE,
  rev(colors_us),
  is_viridis = FALSE,
  legend_off = TRUE
)

enreopy_us = plots_features_compare_us$plot_entropy_0.25 + ylim(-1, 0.1)
peaks_us = plots_features_compare_us$plot_peaks + ylim(-20, 3)
size_us = plots_features_compare_us$plot_int.size + ylim(-3.5, 0.3)

# normalize (normalize = TRUE)
ggsave(paste0(current_path, "entropy_us.pdf"), 
       enreopy_us, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "peaks_us.pdf"), 
       peaks_us,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "int_size_us.pdf"), 
       size_us,
       height = 10, width = 10, units = "cm", dpi=300)

# do not normalize (normalize = FALSE)
error_us = plots_features_compare_us$plot_error + ylim(0.5, 2.3)

ggsave(paste0(current_path, "error_us.pdf"), 
       error_us,
       height = 10, width = 10, units = "cm", dpi=300)


# india
plots_features_compare_in = plot_features_2datasets(
  features_compare_india, 
  vars_to_sum,
  normalize = TRUE,
  rev(colors_in),
  is_viridis = FALSE,
  legend_off = TRUE
)

# normalize (normalize = TRUE)
enreopy_in = plots_features_compare_in$plot_entropy_0.25 + ylim(-1, 0.1)
peaks_in = plots_features_compare_in$plot_peaks + ylim(-20, 3)
size_in = plots_features_compare_in$plot_int.size + ylim(-3.5, 0.3)


ggsave(paste0(current_path, "entropy_india.pdf"), 
       enreopy_in, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "peaks_india.pdf"), 
       peaks_in,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "int_size_india.pdf"), 
       size_in,
       height = 10, width = 10, units = "cm", dpi=300)

# do not normalize (normalize = FALSE)
error_in = plots_features_compare_in$plot_error + ylim(0.5, 2.3)
ggsave(paste0(current_path, "error_india.pdf"), 
       error_in,
       height = 10, width = 10, units = "cm", dpi=300)


# supplementary analysis
colors_culture = c("#54C568FF", "grey")


# ACROSS
plots_features_compare_btw = plot_features_2datasets(
  features_compare_btw, 
  vars_to_sum,
  normalize = TRUE,
  colors_culture, 
  is_viridis = FALSE,
  legend_off = TRUE
)

enreopy_btw = plots_features_compare_btw$plot_entropy_0.25 + ylim(-1, 0.1)
peaks_btw = plots_features_compare_btw$plot_peaks + ylim(-20, 3)
size_btw = plots_features_compare_btw$plot_int.size + ylim(-3.5, 0.3)

# normalize (normalize = TRUE)
ggsave(paste0(current_path, "SI_btw_cross_entropy.pdf"), 
       enreopy_btw, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "SI_btw_cross_peaks.pdf"), 
       peaks_btw,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "SI_btw_cross_size.pdf"), 
       size_btw,
       height = 10, width = 10, units = "cm", dpi=300)


# do not normalize (normalize = FALSE)
error_btw = plots_features_compare_btw$plot_error + ylim(0.5, 2.3)

ggsave(paste0(current_path, "SI_btw_cross_error.pdf"), 
       error_btw,
       height = 10, width = 10, units = "cm", dpi=300)


# WITHIN
plots_features_compare_wth = plot_features_2datasets(
  features_compare_wth, 
  vars_to_sum,
  normalize = FALSE,
  colors_culture, 
  is_viridis = FALSE,
  legend_off = TRUE
)


# normalize (normalize = TRUE)
enreopy_wth = plots_features_compare_wth$plot_entropy_0.25 + ylim(-1, 0.1)
peaks_wth = plots_features_compare_wth$plot_peaks + ylim(-20, 3)
size_wth = plots_features_compare_wth$plot_int.size + ylim(-3.5, 0.3)

ggsave(paste0(current_path, "SI_wth_cross_entropy.pdf"), 
       enreopy_wth, 
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "SI_wth_cross_peaks.pdf"), 
       peaks_wth,
       height = 10, width = 10, units = "cm", dpi=300)
ggsave(paste0(current_path, "SI_wth_cross_size.pdf"), 
       size_wth,
       height = 10, width = 10, units = "cm", dpi=300)


# do not normalize (normalize = FALSE)
error_wth = plots_features_compare_wth$plot_error + ylim(0.5, 2.3)

ggsave(paste0(current_path, "SI_wth_cross_error.pdf"), 
       error_wth,
       height = 10, width = 10, units = "cm", dpi=300)


################################################################################
# stats
################################################################################
boot_stats_exp11 = run_lm_bootraped_data(features_exp11)
boot_stats_exp12 = run_lm_bootraped_data(features_exp12)

# stats 
boot_stats_exp11[[2]]
boot_stats_exp12[[2]]

################################################################################
# Supplementary
## musicians vs non musicians (Experiment 11)
MAX_INTERVAL_SIZE = 13
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

length(table(data_exp11_us.wth$participant_id))  # 184

MT_cut = 3  # 3 years of formal musical training
data_exp11_musicians_wide = data_exp11_us.wth %>% 
  filter(mt > MT_cut)
data_exp11_non.musicians_side = data_exp11_us.wth %>% 
  filter(mt <= MT_cut)

length(table(data_exp11_musicians_wide$participant_id))  # 43
length(table(data_exp11_non.musicians_side$participant_id))  # 141


data_exp11_musicians = data_exp11_musicians_wide %>% 
  pivot_longer(sung_interval1:sung_interval2, values_to = "interval")

data_exp11_musicians %>% 
  distinct(participant_id, .keep_all = T) %>% 
  summarise(mean = mean(mt), sd = sd(mt))

data_exp11_non.musicians = data_exp11_non.musicians_side %>% 
  pivot_longer(sung_interval1:sung_interval2, values_to = "interval")

data_exp11_non.musicians %>% 
  distinct(participant_id, .keep_all = T) %>% 
  summarise(mean = mean(mt), sd = sd(mt))

NBOOT = 1000
marginals_exp11_musicians = make_marginals_kde(data_exp11_musicians, c("interval"), NBOOT, BW)
marginals_exp11_non.musicians = make_marginals_kde(data_exp11_non.musicians, c("interval"), NBOOT, BW)

plots_compare_musicality = plot_grid(marginals_exp11_musicians, marginals_exp11_non.musicians, ncol = 1)

ggsave(paste0(current_path, "SI_marginals_musicality.pdf"), 
       plots_compare_musicality, 
       height = 10, width = 20,
       units = "cm",
       dpi=300)
