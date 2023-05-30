################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis long-melodies (Experiment 3-4)
################################################################################

# Structure: (1) features, (2) contours + clustering, (3) free notes, (4) stats
# Experiments 2-3


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
source("utils/features.R")  

# global parameters
set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

MAX_INTERVAL_SIZE = 12
interval_range = c(-MAX_INTERVAL_SIZE,MAX_INTERVAL_SIZE)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

current_path = "analysis/2-long-melodies/"

data_exp2 <- read_csv("data/experiment2/exp2_btw_5tones_159chains.csv")
data_exp3 <- read_csv("data/experiment3/exp3_btw_free.notes_216chains.csv")


################################################################################
# 1. Features
################################################################################
# peaks from matlab procedure
path_peaksexp1 = "peaks/number_peaks_iterations-btw-3tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp1 = read_boot_peaks_data(path_peaksexp1) 

path_peaksexp2 = "peaks/number_peaks_iterations-btw-5tones-all-intervals.NBOOT-1000.txt"
boot_peaks_exp2 = read_boot_peaks_data(path_peaksexp2) 

path_peaksexp3 = "peaks/number_peaks_iterations-btw-free-all-intervals.NBOOT-1000.txt"
boot_peaks_exp3 = read_boot_peaks_data(path_peaksexp3) 


# features
features_exp1 <- read_csv("features/exp1_features_btw_3tones_590chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp1)

features_exp2 <- read_csv("features/exp2_features_btw_5tones_160chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp2)

features_exp3 <- read_csv("features/exp3_features_btw_free.tones_216chians_1Kboot.csv") %>% 
  mutate(degree=as.factor(degree)) %>% 
  # add peaks
  bind_cols(boot_peaks_exp3)

features_exp1$app = "3-notes"
features_exp2$app = "5-notes"
features_exp3$app = "free-notes"

features_compare = as_tibble(rbind(features_exp1, features_exp2, features_exp3))

# colors
library(scales)
show_col(viridis_pal(option = "C")(2))

# 3 groups
ccc = c("grey", "#0D0887FF", "#F0F921FF")


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
  ccc, 
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


################################################################################
# 2. Contours (Experiment 2)
################################################################################
list_4int_pitches = c("sung_pitch1","sung_pitch2","sung_pitch3", "sung_pitch4", "sung_pitch5")

data_exp2_contour = prepare_data_av.contours(data_exp2, list_4int_pitches)
plot_exp2_contour = plot_av.contour(data_exp2_contour)

ggsave(paste0(current_path, "contours-5tones.pdf"), 
       plot_exp2_contour,
       height = 8, width = 14, units = "cm", dpi=300)


# proportion intervals larger than 7
list_4int_sung.intervals = c("sung_interval1","sung_interval2","sung_interval3","sung_interval4")

proportion_large_ints = get_proportion_large_ints_boot(data_exp2, list_4int_sung.intervals, 10)

options(pillar.sigfig = 4)
proportion_large_ints %>% 
  group_by(degree) %>% 
  dplyr::summarise(
    n = n(),
    mean = mean(proportion),
    sd = sd(proportion) + 1.96,  # to get 95 CI
    low_ci = mean - sd,
    high_ci = mean + sd
  )

# clustering
list_4int_pitches = c("sung_pitch1","sung_pitch2","sung_pitch3", "sung_pitch4", "sung_pitch5")
list_4int_int2refs =  c("in2ref_1", "in2ref_2", "in2ref_3", "in2ref_4", "in2ref_5")

data_to.clust_exp2 = data_exp2 %>%  
  filter(degree %in% 8:10) %>% 
  select(network_id, sung_pitch1:sung_pitch5) %>% 
  mutate(mean_pitch = (sung_pitch1+sung_pitch2+sung_pitch3+sung_pitch4+sung_pitch5)/5) %>%
  mutate(
    in2ref_1 = sung_pitch1 - mean_pitch,
    in2ref_2 = sung_pitch2 - mean_pitch,
    in2ref_3 = sung_pitch3 - mean_pitch,
    in2ref_4 = sung_pitch4 - mean_pitch,
    in2ref_5 = sung_pitch5 - mean_pitch
  ) %>% 
  select(-sung_pitch1:-sung_pitch5, -mean_pitch) 


# clustering
library(factoextra)
library(cluster)
library(ggfortify)
library(stats)

fviz_nbclust(data_to.clust_exp2[,-1], kmeans, method = "wss")
k = 2

pca_res = prcomp(data_to.clust_exp2[,-1])
summary(pca_res)
cluster <- kmeans(data_to.clust_exp2[,-1], k, nstart = 50)

data_to.clust_exp2$cluster <- as.factor(cluster$cluster)


show_col(viridis_pal(option = "F")(15))
ccc2 = c("#2A788EFF", "#E13342FF")

pca_plot = autoplot(pca_res,  data = data_to.clust_exp2, colour = 'cluster', size = 1) +
  scale_color_manual(values=ccc2)  +
  theme(axis.text.x = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=10),
        legend.title = element_blank())

plot_cluster_solutions = plot_clustered_contours(data_to.clust_exp2, k, list_4int_int2refs, ccc2) 


ggsave(paste0(current_path, "pca_plot.pdf"), 
       pca_plot,
       height = 10, width = 10, units = "cm", dpi=300)

ggsave(paste0(current_path, "plot_cluster_solutions.pdf"), 
       plot_cluster_solutions,
       height = 8, width = 8, units = "cm", dpi=300)


################################################################################
# 3. Free notes (Experiment 3)
################################################################################
plot_overall_trend = data_exp3 %>% 
  # plot_free_n.notes_overall.trend("#CC4678FF")
  plot_free_n.notes_overall.trend("#F0F921FF")

plot_final_hist = data_exp3 %>% 
  filter(num_sung_pitches %in% 4:12) %>% 
  plot_free_n.notes_last.histogram()

ggsave(paste0(current_path, "free.notes_trend.pdf"), 
       plot_overall_trend,
       units = "cm",
       width = 8, height = 8, dpi = 300) 
ggsave(paste0(current_path, "free.notes_histogram.pdf"), 
       plot_final_hist,
       units = "cm",
       width = 7, height = 7, dpi = 300) 


# contours
# contours
data_contours = data_exp3 %>% 
  mutate(num_sung_pitches = ifelse(degree == 0, num_target_pitches, num_sung_pitches)) %>% 
  mutate(sung_pitches = map(strsplit(sung_pitches, ","), parse_number)) %>% 
  filter(degree %in% 8:10)

cols = viridis_pal(option = "D")(8)

plot_contours = plot_grid(
  plot_av_contours_general(data_contours, 4, cols[2]) + ggtitle("4 notes"),
  plot_av_contours_general(data_contours, 5, cols[3]) + ggtitle("5 notes"),
  plot_av_contours_general(data_contours, 6, cols[4]) + ggtitle("6 notes"),
  plot_av_contours_general(data_contours, 7, cols[5]) + ggtitle("7 notes"),
  plot_av_contours_general(data_contours, 8, cols[6]) + ggtitle("8 notes"),
  plot_av_contours_general(data_contours, 9, cols[7]) + ggtitle("9 notes")
)


ggsave(paste0(current_path, "free.notes_contours.pdf"), 
       plot_contours, 
       units = "cm",
       width = 17, height = 12, dpi = 300) 


################################################################################
# 4. Free notes (Experiment 3)
################################################################################

boot_stats_exp2 = run_lm_bootraped_data(features_exp2)
boot_stats_exp3 = run_lm_bootraped_data(features_exp3)

# stats 
boot_stats_exp2[[2]]
boot_stats_exp3[[2]]


