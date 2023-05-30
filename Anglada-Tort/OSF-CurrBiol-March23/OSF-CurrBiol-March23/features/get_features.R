################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Melodic features via bootstrap
################################################################################
# imports
library(tidyverse)

# methods
source("utils/features.R")  # methods for calculating melodic features

# data melodic consonance
data_melcon_smoothed <- read_csv("data/experiment6/exp6_melodic.pleasant_all.batches_smoothed.csv")

# Global parameters
interval_range = c(-15,15)
vertical.lines = seq(from=min(interval_range), to=max(interval_range), by = 1)

list_2int_sung.intervals = c("sung_interval1","sung_interval2")
BW = 0.25
NBOOT = 1000

current_path = "features/"


################################################################################
# Experiment 1: 3-note melodies (across)
################################################################################
# data Experiment 1 (across)
exp1_btw_3tones_590chains <- read_csv("data/experiment1/exp1_btw_3tones_590chains.csv")

length(unique(exp1_btw_3tones_590chains$participant_id)) # 188
length(unique(exp1_btw_3tones_590chains$network_id)) # 600
table(exp1_btw_3tones_590chains$degree) # 590

features_exp1_btw_3tones_590chains = get_bootstrapped_features(
  exp1_btw_3tones_590chains, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT
  )

name_file = "exp1_features_btw_3tones_590chians_1Kboot"

features_exp1_btw_3tones_590chains %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 2: 5-note melodies
################################################################################
exp2_btw_5tones_all <- read_csv("data/experiment2/exp2_btw_5tones_159chains.csv")
list_4int_sung.intervals = c("sung_interval1","sung_interval2","sung_interval3","sung_interval4")

length(unique(exp2_btw_5tones_all$participant_id)) # 51
length(unique(exp2_btw_5tones_all$network_id)) # 160
table(exp2_btw_5tones_all$degree) # 159

exp2_features_btw_5tones_all = get_bootstrapped_features(
  exp2_btw_5tones_all, 
  data_melcon_smoothed,
  list_4int_sung.intervals,
  NBOOT
)

name_file = "exp2_features_btw_5tones_160chians_1Kboot"

exp2_features_btw_5tones_all %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 3: 4-12 note melodies
################################################################################
exp3_btw_free_all <- read_csv("data/experiment3/exp3_btw_free.notes_216chains.csv") %>% 
  mutate(sung_intervals = map(strsplit(sung_intervals, ","), parse_number)) %>% 
  rename(list_sung_ints=sung_intervals)

length(unique(exp3_btw_free_all$participant_id)) # 83
length(unique(exp3_btw_free_all$network_id)) # 216
table(exp3_btw_free_all$degree) # 216


exp3_features_btw_free_all = get_bootstrapped_features(
  exp3_btw_free_all, 
  data_melcon_smoothed,
  list_4int_sung.intervals,
  NBOOT,
  is_free_notes = TRUE
)

name_file = "exp3_features_btw_free.tones_216chians_1Kboot"

exp3_features_btw_free_all %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 4 & 5: 2 notes singing vs slider
################################################################################
# slider
exp4_slider <- read_csv("data/experiment4/exp4_btw_2tones_slider_369chains.csv") %>% 
  distinct(origin_id, .keep_all = TRUE) %>%   # remove duplicated from aggregation
  dplyr::rename(interval=location)

length(unique(exp4_slider$participant_id)) # 327
length(unique(exp4_slider$network_id)) # 400
table(exp4_slider$degree) # 369

exp4_features_btw_slider_2tones_all = get_features_bootstraped_slider(
  exp4_slider, 
  data_melcon_smoothed,
  c("interval"),
  NBOOT
)

name_file = "exp4_features_btw_slider_2tones_372chains_1Kboot"

exp4_features_btw_slider_2tones_all %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


# singing
exp5_2tones <- read_csv("data/experiment5/exp5_btw_2tones_398chains.csv")  %>% 
  dplyr::rename(interval = sung_interval1) %>% 
  mutate(interval_error = abs(interval - target_interval1))

length(unique(exp5_2tones$participant_id)) # 122
length(unique(exp5_2tones$network_id)) # 400
table(exp5_2tones$degree) # 398

exp5_features_btw_2tones_all = get_features_bootstraped_slider(
  exp5_2tones, 
  data_melcon_smoothed,
  c("interval"),
  NBOOT
)

name_file = "exp5_features_btw_2tones_398chains_1Kboot"

exp5_features_btw_2tones_all %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 7-9: memory 
################################################################################
data_exp7 <- read_csv("data/experiment7-9/exp7_btw_3tones_memory.5sec_240chains.csv")
data_exp8 <- read_csv("data/experiment7-9/exp8_btw_3tones_memory.10sec_240chains.csv")
data_exp9 <- read_csv("data/experiment7-9/exp9_btw_3tones_memory.control_240chains.csv")

# memory 5 sec
length(unique(data_exp7$participant_id)) # 100
length(unique(data_exp7$network_id)) # 240
table(data_exp7$degree) # 240


features_data_exp7 = get_bootstrapped_features(
  data_exp7, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT, 
)

name_file = "exp7_features_btw_3tones_memory_240chians_1Kboot"

features_data_exp7 %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


# memory 10 sec
length(unique(data_exp8$participant_id)) # 105
length(unique(data_exp8$network_id)) # 240
table(data_exp8$degree) # 240


features_data_exp8 = get_bootstrapped_features(
  data_exp8, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT, 
)

name_file = "exp8_features_btw_3tones_memory_240chians_1Kboot"

features_data_exp8 %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


# memory (control)
length(unique(data_exp9$participant_id)) # 95
length(unique(data_exp9$network_id)) # 240
table(data_exp9$degree) # 240


features_data_exp9 = get_bootstrapped_features(
  data_exp9, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT, 
)

name_file = "exp9_features_btw_3tones_memory_240chians_1Kboot"

features_data_exp9 %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 10: INDIA, 3-note melodies (across)
################################################################################
exp10_btw_3tones_india_120chains <- read_csv("data/experiment10/exp10_btw_3tones_india_120chains.csv") %>% 
  group_by(network_id) %>% 
  mutate(n_trials = n()) %>% 
  filter(n_trials > 9)  # completing all chains in this experiment was not possible because the limited pool of active participants


length(unique(exp10_btw_3tones_india_120chains$participant_id)) # 54
length(unique(exp10_btw_3tones_india_120chains$network_id)) # 200
table(exp10_btw_3tones_india_120chains$degree) # 120 (9th generation)

features_exp10_btw_3tones_india_120chains = get_bootstrapped_features(
  exp10_btw_3tones_india_120chains, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT
)

name_file = "exp10_features_btw_3tones_india_120chains_1Kboot"

features_exp10_btw_3tones_india_120chains %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


################################################################################
# Experiment 11: US, 3-note melodies (within)
################################################################################
exp11_wth.us_3tones_all <- read_csv("data/experiment11/exp11_wth_3tones_615chains.csv")

length(unique(exp11_wth.us_3tones_all$participant_id)) # 162
length(unique(exp11_wth.us_3tones_all$network_id)) # 615
table(exp11_wth.us_3tones_all$degree) # 615

features_exp11_wth.us_3tones_all = get_bootstrapped_features(
  exp11_wth_3tones_all, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT
)

name_file = "exp11_features_wth_3tones_615chians_1Kboot"

exp11_features_wth_3tones_all %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)



################################################################################
# Experiment 12: 3 notes INDIA (within) 
################################################################################
exp12_wth_3tones_india_223chains <- read_csv("data/experiment12/exp12_wth_3tones_india_223chains.csv") %>% 
  group_by(network_id) %>% 
  mutate(n_trials = n()) %>% 
  filter(n_trials > 10)  # only accept full chains


length(unique(exp12_wth_3tones_india_223chains$participant_id)) # 74
length(unique(exp12_wth_3tones_india_223chains$network_id)) # 223
table(exp12_wth_3tones_india_223chains$degree) # 223


features_exp12_wth_3tones_india_223chains = get_bootstrapped_features(
  exp12_wth_3tones_india_223chains, 
  data_melcon_smoothed,
  list_2int_sung.intervals, 
  NBOOT
)

name_file = "exp12_features_wth_3tones_india_223chains_1Kboot"

features_exp12_wth_3tones_india_223chains %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv(paste0(current_path, name_file,".csv"), row.names = FALSE)


