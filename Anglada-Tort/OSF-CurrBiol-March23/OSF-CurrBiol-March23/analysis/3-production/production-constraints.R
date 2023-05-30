################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis production constraints (Experiment 5-6)
################################################################################

# Structure: (1) distributions, (2) features, (3) preferences, (4) stats
# Experiment 4-6

################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(egg)
library(ggpubr)
library(cowplot)

# methods
source("utils/plots.R")  # methods for plotting
source("utils/stats.R")  # general methods

# global parameters
set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

MAX_INTERVAL_SIZE = 13
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

current_path = "analysis/3-production/"

# colors
library(scales)
show_col(viridis_pal(option = "A")(10))
#  groups
colors_slider = c("grey", "#E13342FF")


################################################################################
# Distributions
################################################################################
# slider
exp4_slider <- read_csv("data/experiment4/exp4_btw_2tones_slider_369chains.csv") %>% 
  distinct(origin_id, .keep_all = TRUE) %>%   # remove duplicated from aggregation
  dplyr::rename(interval=location)
# singing
exp5_2tones <- read_csv("data/experiment5/exp5_btw_2tones_398chains.csv")  %>% 
  dplyr::rename(interval = sung_interval1)


length(unique(exp4_slider$participant_id)) # 327
table(exp4_slider$degree)  # 369

length(unique(exp5_2tones$participant_id)) # 122
table(exp5_2tones$degree)  # 398


# marginals
NBOOT = 1000
BW = 0.25


# together
peaks_exp4_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-slider-2tones-all-intervals.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

peaks_exp4 = get_sig_peaks(peaks_exp4_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range))


marginals_slider_singing = make_marginals_2datasets_slider(
  exp5_2tones, exp4_slider,
  peaks_exp4,
  c("interval"),
  NBOOT, BW, interval_range,
  colors_slider
)

ggsave(paste0(current_path, "marginals_slider_singing.pdf"), 
       marginals_slider_singing, 
       width = 20, height = 8, units = "cm", dpi=300, bg = "white")


################################################################################
# Supplementary

marginals_2tones_slider = make_marginals_kde(exp4_slider, c("interval"), NBOOT, BW)
marginals_2tones = make_marginals_kde(exp5_2tones, c("interval"), NBOOT, BW)


ggsave(paste0(current_path, "SI_marginals_2tones_slider.pdf"), 
       marginals_2tones_slider, 
       height = 10, width = 20,
       units = "cm",
       dpi=300)


ggsave(paste0(current_path, "SI_sing_marginals_2tones.pdf"), 
       marginals_2tones, 
       height = 10, width = 20,
       units = "cm",
       dpi=300)


################################################################################
# Features
################################################################################
# peaks from matlab procedure
path_peaks = "peaks/number_peaks_iterations-slider-2tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp4 = read_boot_peaks_data(path_peaks) 

path_peaks = "peaks/number_peaks_iterations-btw-2tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp5 = read_boot_peaks_data(path_peaks) 

# features
exp4_features_2tones_slider <- read_csv("features/exp4_features_btw_slider_2tones_372chains_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp4)
exp5_features_2tones <- read_csv("features/exp5_features_btw_2tones_398chains_1Kboot.csv") %>% 
  mutate(degree = as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp5)

exp4_features_2tones_slider$app = "slider-2t"
exp5_features_2tones$app = "sing-2t"


features_compare = as_tibble(rbind(exp4_features_2tones_slider, exp5_features_2tones)) %>% 
  mutate(mean_error_interval = ifelse(degree == 0, NA, mean_error_interval))

vars_to_sum = c(
  "mean_error_interval",
  "mean_abs_int_size",
  "interval_entropy_0.25",
  "n_peaks"
)

plots_features_compare = plot_features_2datasets_slider(
  features_compare, 
  vars_to_sum,
  normalize = TRUE,
  colors_slider, 
  is_viridis = FALSE,
  legend_off = TRUE
)

plots_features_compare$plot_entropy_0.25
plots_features_compare$plot_peaks
plots_features_compare$plot_int.size

plots_features_compare$plot_error


# normalize (normalize = TRUE)
ggsave(paste0(current_path, "entropy_slider.pdf"), 
       plots_features_compare$plot_entropy_0.25, 
       width = 10, height = 10, units = "cm", dpi=300, bg = "white")
ggsave(paste0(current_path, "peaks_slider.pdf"), 
       plots_features_compare$plot_peaks, 
       width = 10, height = 10, units = "cm", dpi=300, bg = "white")
ggsave(paste0(current_path, "int_size_slider.pdf"), 
       plots_features_compare$plot_int.size, 
       width = 10, height = 10, units = "cm", dpi=300, bg = "white")

# do not normalize (normalize = FALSE)
ggsave(paste0(current_path, "error.pdf"), 
       plots_features_compare$plot_error, 
       width = 10, height = 10, units = "cm", dpi=300, bg = "white")


################################################################################
# Melodic pleasantness (Experiment 6)
################################################################################
interval_range = c(-15, 15)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

# rating experiment
exp6_data_melcon_smoothed <- read_csv("data/experiment6/exp6_melodic.pleasant_all.batches_smoothed.csv")


# ## Comment out to smoth raw ratings (it is slow)
# exp6_data_melcon_raw <- read_csv("data/experiment6/exp6_melodic.pleasant_all.batches.csv")
# # function to smooth ratings
# Rcpp::sourceCpp("currentbio22/paper/utils/smooth_1d_gaussian.cpp")
# # smooth ratings (1K boots)
# exp6_data_melcon_smoothed = smooth_melcon_ratings(exp6_data_melcon_raw)
# 
# # save smoothed data
# write_csv(exp6_data_melcon_smoothed, "currentbio22/paper/data/experiment6/exp6_melodic.pleasant_all.batches_smoothed.csv")
# 
# # END 
# ##


# plot melodic pleasantness
interval_range = c(-13, 13)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

plot_rating_continous = exp6_data_melcon_smoothed %>% 
  filter(interval <= 13 & interval >= -13) %>% 
  plot_melcon_ratings("#35B779FF", interval_range) 


# add peaks to plot
peaks_exp6_all = as_tibble(read.table(
  "peaks/boottsrapping_peaks_matrix-melcon-2tones-all.NBOOT-1000.txt",
  header = FALSE, sep = ",", dec = "."))

peaks_exp6 = get_sig_peaks(peaks_exp6_all, 90) %>%  
  filter(values > min(interval_range) & values < max(interval_range)) %>% 
  filter(values <= 13 & values >= -13) 

# get y coordinates based on x
plot_rating_continous_ys = plot_melcon_ratings_onlyline(exp6_data_melcon_smoothed, interval_range)
peaks_exp6_x.y = get_sing_dots_marginals_ratings(plot_rating_continous_ys, peaks_exp6)


plot_rating_continous_peaks = add_peaks_marginals_loop(
  plot_rating_continous, 
  peaks_exp6_x.y, 
  peaks_exp6, 
  "red")


ggsave(paste0(current_path, "melodic_pleasant_profile_with.singing.pdf"), 
       plot_rating_continous_peaks, 
       height = 8, width = 20,
       units = "cm",
       dpi=300)

################################################################################
# Singing alignment with melodic pleasantness
################################################################################
features_exp1 <- read_csv("features/exp1_features_btw_3tones_590chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree)) %>% 
  select(boot, degree, melcon)
features_exp2 <- read_csv("features/exp2_features_btw_5tones_160chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree))  %>% 
  select(boot, degree, melcon)
features_exp3 <- read_csv("features/exp3_features_btw_free.tones_216chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree)) %>% 
  select(boot, degree, melcon)

exp5_features_2tones_melcon =  exp5_features_2tones %>% 
  select(boot, degree, melcon)

exp5_features_2tones_melcon$app = "2t"
features_exp1$app = "3t"
features_exp2$app = "5t"
features_exp3$app = "freet"

features_compare = as_tibble(rbind(
  exp5_features_2tones_melcon,
  features_exp1,features_exp2, features_exp3)) 

features_compare_sum_pre = sum_boot.data_by.group(features_compare, "melcon") %>% 
  group_by(app) %>% 
  mutate(melcon = mean-mean[degree==0]) # normalize by first generation

show_col(viridis_pal(option = "D")(4))
colors_melcon = c("grey", "#440154FF", "#31688EFF", "#FDE725FF")

plot_melcon = make_lineplot_ribon_bygroup(
  features_compare_sum_pre, melcon, sd, app, 
  "Melodic Pleasantness",
  colors_melcon, FALSE, FALSE) 

ggsave(paste0(current_path, "plot_melcon.pdf"), 
       plot_melcon, 
       height = 10, width = 12,
       units = "cm",
       dpi=300)


################################################################################
# stats
################################################################################
exp4_features_2tones_slider_renamed = exp4_features_2tones_slider %>% 
  # renaming features as required for the input names in the stats funcion
  rename(RMSE_interval = mean_error_interval, abs_int_size = mean_abs_int_size,
         entropy_0.25 = interval_entropy_0.25)

exp5_features_2tones_renamed = exp5_features_2tones %>% 
  # renaming features as required for the input names in the stats funcion
  rename(RMSE_interval = mean_error_interval, abs_int_size = mean_abs_int_size,
         entropy_0.25 = interval_entropy_0.25) %>% 
  mutate(RMSE_interval = ifelse(degree==0, NA, RMSE_interval))

boot_stats_exp4 = run_lm_bootraped_data(exp4_features_2tones_slider_renamed)
boot_stats_exp5 = run_lm_bootraped_data(exp5_features_2tones_renamed)

# stats 
boot_stats_exp4[[2]]
boot_stats_exp5[[2]]

