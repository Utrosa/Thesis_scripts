################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Power Analysis
################################################################################

# Structure: (1) sample-size analysis Experiment 1 and (2) Experiment 6

################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(performance)  # splithalf correlation
library(philentropy)  # JSD
library(mgcv)  # gam
library(cowplot)

current_path = "sample-size/"

source("utils/plots.R")  # methods for plotting
source("utils/stats.R")  # methods for stats
source("utils/features.R")  # methods for calculating melodic features

#global params
set.seed(2022)
BW = 0.25
interval_range = c(-15,15)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

################################################################################
# run power analysis Experiment 1
################################################################################
N_list = seq(50, 600, 50)

# Set FALSE to import data; set TRUE to run the power analysis again (takes long)
RUN_POWER_ANALYSIS = FALSE


if (RUN_POWER_ANALYSIS){
  print("run power analysis (takes a long time)")
  # data
  data_exp1 <- read_csv("data/experiment1/exp1_btw_3tones_590chains.csv")
  
  res_sample.size_exp1 = run_sample.size_analysis(
    data_exp1, 
    BW,
    N_list,  
    N_boot = 1000  # this should be 1K
  )
  
  res_sample.size_exp1
  
  # data to plot
  results_sample.size_sum_exp1 = summarize_results_boot(res_sample.size_exp1) %>% 
    filter(metric == "JSD" | metric == "peaks_spearmanbrown" | metric == "entropy.cor.r")
  
  # save
  res_sample.size_exp1.save = res_sample.size_exp1 %>%  select(-output)
  write_csv(res_sample.size_exp1.save, "results_sample-size.analysis_exp1_1Kboots.csv")
  
} else {
  print("import data")
  res_sample.size_exp1 <- read_csv("sample-size/results_sample-size.analysis_exp1_1Kboots.csv")
  
  # data to plot
  results_sample.size_sum_exp1 = res_sample.size_exp1 %>% 
    group_by(N) %>% 
    dplyr::summarise(
      n = n(),
      # mean
      JSD=mean(jsd),
      entropy.cor.r=mean(entropy.cor.r),
      entropy.cor.sp=mean(entropy.cor.sp),
      peaks_splithalf=mean(peaks_splithalf),
      peaks_spearmanbrown=mean(peaks_spearmanbrown)
    ) %>% 
    select(N, JSD:peaks_spearmanbrown) %>% 
    pivot_longer(JSD:peaks_spearmanbrown, "metric", values_to = "mean")  %>% 
    filter(metric == "JSD" | metric == "peaks_spearmanbrown" | metric == "entropy.cor.r")
  
}


# default colors in ggplot
library(scales)
show_col(hue_pal()(3))


plot_jsd = plot_power.analysis_metric(
  results_sample.size_sum_exp1, "JSD", "#f8766D", "JSD (splithalf)",
  "Jensen-Shannon divergence"
) +
  scale_x_continuous(breaks=seq(min(N_list), max(N_list), 50)) +
  scale_y_continuous(breaks=seq(0, 0.15, 0.02)) +
  # add line where distribution becomes flat
  geom_hline(aes(yintercept = 0.01), linetype="dashed", color = "red") 


plot_entropy = plot_power.analysis_metric(
  results_sample.size_sum_exp1, "entropy.cor.r", "#619CFF", "r (splithalf)",
  "Correlation peaks (rkk)"
) + 
  # add line where r is high (r = 0.8)
  geom_hline(aes(yintercept = 0.8), linetype="dashed", color = "red") +
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0,1)) +
  scale_x_continuous(breaks=seq(min(N_list), max(N_list), 50)) 

power_analysis_plots = plot_grid(plot_jsd, plot_entropy, nrow = 1)


ggsave(paste0(current_path, "power_analysis_plots_exp1.pdf"), 
       power_analysis_plots, 
       width = 20, height = 10, units = "cm", dpi=300, bg = "white")


################################################################################
# power analysis Experiment 6
################################################################################
N_list_ratings = seq(500, 15000, 500)

# Set FALSE to import data; set TRUE to run the power analysis again (takes long)
RUN_POWER_ANALYSIS = FALSE


if (RUN_POWER_ANALYSIS){
  print("run power analysis (takes a long time)")
  # data
  
  Rcpp::sourceCpp("utils/smooth_1d_gaussian.cpp")
  
  # data
  data_exp6 <- read_csv("data/experiment6/exp6_melodic.pleasant_all.batches.csv") %>% 
    filter(is_repeat_trial == FALSE) %>% 
    select(id, participant_id, network_id, stimulus_id, answer, target_interval)
  
  length(table(data_exp6$stimulus_id))
  
  res_sample.size_exp6 = run_sample.size_analysis_ratings(
    data_exp6, 
    BW,
    N_list_ratings,  
    N_boot = 1000  # this should be 1K
  )
  
  res_sample.size_exp6
  
  # save
  res_sample.size_exp6.save = res_sample.size_exp6 %>%  select(-output)
  write_csv(res_sample.size_exp6.save, "results_sample.size.analysis_exp6_1Kboots.csv")
  
  # data to plot
  results_sample.size_sum_exp6 = summarize_results_boot(res_sample.size_exp6) %>% 
    filter(metric == "JSD" | metric == "peaks_spearmanbrown") 
  
} else {
  print("import data")
  res_sample.size_exp6 <- read_csv("sample-size/results_sample-size.analysis_exp6_1Kboots.csv")
  
  # data to plot
  results_sample.size_sum_exp6 = res_sample.size_exp6 %>% 
    group_by(N) %>% 
    dplyr::summarise(
      n = n(),
      # mean
      JSD=mean(jsd),
      peaks_splithalf=mean(peaks_splithalf),
      peaks_spearmanbrown=mean(peaks_spearmanbrown)
    ) %>% 
    select(N, JSD:peaks_spearmanbrown) %>% 
    pivot_longer(JSD:peaks_spearmanbrown, "metric", values_to = "mean")  %>% 
    filter(metric == "JSD" | metric == "peaks_spearmanbrown" | metric == "entropy.cor.r")
  
}


plot_jsd.ratings = plot_power.analysis_metric(
  results_sample.size_sum_exp6, "JSD",  "#00BA38", "JSD (splithalf)",
  "Jensen-Shannon divergence"
) +
  scale_x_continuous(breaks=seq(min(N_list_ratings), max(N_list_ratings), 2000)) +
  scale_y_continuous(breaks=seq(0, 0.02, 0.001)) +
  # add line where distribution becomes flat
  geom_hline(aes(yintercept = 0.0003), linetype="dashed", color = "red")


ggsave(paste0(current_path, "power_analysis__exp6.pdf"), 
       plot_jsd.ratings, 
       width = 10, height = 10, units = "cm", dpi=300, bg = "white")

