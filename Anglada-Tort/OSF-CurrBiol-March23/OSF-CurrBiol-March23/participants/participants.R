################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Participant demographcis and singing accuracy
################################################################################

# Structure: (1) basic demographics, (2) musical experience, (3) singing accuracy


################################################################################
# Global
################################################################################
# imports
library(tidyverse)

current_path = "currentbio22/paper/participants/"


################################################################################
# 1. demographics
################################################################################

# summary functions
get_demographics = function(data){
  unique.pids.data <- unique(data[, c("participant_id", "gender")]) %>% 
    filter(gender == "female" | gender == "male")
  table_gender = table(unique.pids.data$gender)
  demographics = tibble(
    N = length(table(data$participant_id)),
    mean_age = mean(data$age, na.rm = T), 
    sd_age = sd(data$age, na.rm = T), 
    min_age = min(data$age, na.rm = T), 
    max_age = max(data$age, na.rm = T), 
    N_female = table_gender[[1]],
    Percent_female =  ( table_gender[[1]]/N)* 100,
    N_male =  table_gender[[2]],
    Percent_male =  ( table_gender[[2]]/N)* 100,
    M_mt = mean(data$mt, na.rm = T), 
    SD_mt = sd(data$mt, na.rm = T), 
    M_GMS_MT = mean(data$GMSI_MT, na.rm = T), 
    SD_GMS_MT = sd(data$GMSI_MT, na.rm = T),
    M_GMS_SA = mean(data$GMSI_SA, na.rm = T), 
    SD_GMS_SA = sd(data$GMSI_SA, na.rm = T)
  )
  return(demographics)
}


calculate_from_raw_batch  = function(data, batch_id){
  res = data %>% 
    dplyr::filter(failed == FALSE) %>% 
    select(id, participant_id, question, response) %>% 
    dplyr::filter( question == "gender" | question == "age" | question == "years_of_formal_training") %>%
    pivot_wider(id_cols = participant_id, names_from = question, values_from = response) %>%
    mutate(
      years_of_formal_training = parse_number(as.character(years_of_formal_training)), 
      age = parse_number(as.character(age)), 
      gender = gsub("[^[:alpha:]]", "", gender)
    ) %>%  
    dplyr::rename(mt=years_of_formal_training) %>% 
    mutate(participant_id = paste0(participant_id, batch_id))
  
  return(res)
}


################################################################################
# Demographic information per experiment
################################################################################
# exp1
data_exp1 <- read_csv("data/experiment1/exp1_btw_3tones_590chains.csv")
table(data_exp1$degree)  # 600/ 590
length(unique(data_exp1$participant_id))  # 188
get_demographics(data_exp1)

# exp2
data_exp2 <- read_csv("data/experiment2/exp2_btw_5tones_159chains.csv")
table(data_exp2$degree)  # 160/ 159
length(unique(data_exp2$participant_id))  # 51
get_demographics(data_exp2)

# exp3
data_exp3 <- read_csv("data/experiment3/exp3_btw_free.notes_216chains.csv")
table(data_exp3$degree)  # 216/ 216
length(unique(data_exp3$participant_id))  # 83
get_demographics(data_exp3)

# exp4
data_exp4 <- read_csv("data/experiment4/exp4_btw_2tones_slider_369chains.csv") %>% 
  distinct(origin_id, .keep_all = TRUE)  # remove aggregated response (only include median)
table(data_exp4$degree)  # 400/ 369
length(unique(data_exp4$participant_id))  # 327
get_demographics(data_exp4)

# exp5
data_exp5 <- read_csv("data/experiment5/exp5_btw_2tones_398chains.csv") 
table(data_exp5$degree)  # 400/ 394
length(unique(data_exp5$participant_id))  # 120
get_demographics(data_exp5)

# exp6
# !! Data was not saved in the main csv, thus using raw psynet export csvs
response_melcon_b1 <- read_csv("data/experiment6/demographics/response_melcon_b1.csv")
response_melcon_b2 <- read_csv("data/experiment6/demographics/response_melcon_b2.csv")
response_melcon_b3 <- read_csv("data/experiment6/demographics/response_melcon_b3.csv")

dem_melcon_b1 = calculate_from_raw_batch(response_melcon_b1, "b1")
dem_melcon_b2 = calculate_from_raw_batch(response_melcon_b2, "b2")
dem_melcon_b3 = calculate_from_raw_batch(response_melcon_b3, "b3")

data_exp6 = as_tibble(rbind(dem_melcon_b1, dem_melcon_b2, dem_melcon_b3)) %>% 
  filter(age < 200)
length(unique(data_exp6$participant_id))  # 416
get_demographics(data_exp6)   # TODO: I need demographics

# exp7-9
data_exp7 <- read_csv("data/experiment7-9/exp7_btw_3tones_memory.5sec_240chains.csv")
table(data_exp7$degree)  # 240/ 240
length(unique(data_exp7$participant_id))  # 105
get_demographics(data_exp7)

data_exp8 <- read_csv("data/experiment7-9/exp8_btw_3tones_memory.10sec_240chains.csv")
table(data_exp8$degree)  # 240/ 240
length(unique(data_exp8$participant_id))  # 100
get_demographics(data_exp8)

data_exp9 <- read_csv("data/experiment7-9/exp9_btw_3tones_memory.control_240chains.csv")
table(data_exp9$degree)  # 240/ 240
length(unique(data_exp9$participant_id))  # 95
get_demographics(data_exp9)

# exp10
data_exp10 <- read_csv("data/experiment10/exp10_btw_3tones_india_120chains.csv")
table(data_exp10$degree)  # 200/ 120-95
length(unique(data_exp10$participant_id))  # 54
get_demographics(data_exp10)

# exp11
data_exp11 <- read_csv("data/experiment11/exp11_wth_3tones_615chains.csv") %>% 
  group_by(network_id) %>%  mutate(n_trials_net = n()) %>% 
  filter(n_trials_net > 9)  # completing this experiment was difficult because the limited pool of active participants
table(data_exp11$degree)  # 615/ 615
length(unique(data_exp11$participant_id))  # 184
get_demographics(data_exp11)

# exp12
data_exp12 <- read_csv("data/experiment12/exp12_wth_3tones_india_223chains.csv")
table(data_exp12$degree)  # 237/ 223
length(unique(data_exp12$participant_id))  # 73
get_demographics(data_exp12)


################################################################################
# 2. singing accuracy (performance test)
################################################################################
get_singing_perforamnce_scores = function(data){
  sing.performance_sum = data %>% 
    group_by(participant_id) %>% 
    dplyr::summarise(
      n = n(),
      mean_int_accuracy = mean(max_abs_interval_error),
      sd_int_accuracy = sd(max_abs_interval_error),
      mean_direction_accuracy = mean(direction_accuracy),
      sd_direction_accuracy = sd(direction_accuracy)
    )
  
  out = tibble(
    min = min(sing.performance_sum$mean_int_accuracy),
    max = max(sing.performance_sum$mean_int_accuracy),
    mean = mean(sing.performance_sum$mean_int_accuracy),
    sd = sd(sing.performance_sum$mean_int_accuracy)
  )
  return(out)
}

prepare_accuracy_data = function(data, batch){
  data_packed = data %>% 
    dplyr::filter(failed == FALSE & type=="singing_performance_trial") %>% 
    select(-creation_time:-time_of_death)
  
  data_packed$analysis[is.na(data_packed$analysis)] <- "{}"
  
  data_unpacked = unpack_json_column(data_packed, data_packed$analysis) %>%
    dplyr::select(id, network_id, participant_id, max_abs_interval_error, direction_accuracy) %>% 
    mutate(participant_id = paste0(participant_id, batch)) %>% 
    mutate(network_id = paste0(network_id, batch)) %>% 
    mutate(id = paste0(id, batch))
  
  return(data_unpacked)
}

sort_json <- function(x){
  jsonlite::stream_in(textConnection(gsub("\\n", "", x)))
}


unpack_json_column = function(data, column_to_unpack){
  column_unpacked = sort_json(column_to_unpack)
  data_unpacked = as_tibble(cbind(data, column_unpacked), .name_repair = "universal")
  return(data_unpacked)
}


# exp1
data_exp1_demographics =
  data_exp1 %>% 
  select(participant_id, age, gender, GMSI_MT, GMSI_SA, mt) %>% 
  distinct(participant_id, .keep_all = T)

exp1_batch1 <- read_csv("data/singing-performance-test/experiment1/exp1_batch1.csv")
exp1_batch2 <- read_csv("data/singing-performance-test/experiment1/exp1_batch2.csv")
exp1_batch3 <- read_csv("data/singing-performance-test/experiment1/exp1_batch3.csv")

exp1_batch1_sing = prepare_accuracy_data(exp1_batch1, "b1")
exp1_batch2_sing = prepare_accuracy_data(exp1_batch2, "b2")
exp1_batch3_sing = prepare_accuracy_data(exp1_batch3, "b3")

exp1_sing.performance = as_tibble(rbind(exp1_batch1_sing, exp1_batch2_sing, exp1_batch3_sing))
unique(exp1_sing.performance[, c("participant_id")])  # 219

get_singing_perforamnce_scores(exp1_sing.performance)


exp1_sing.performance_cor =  exp1_sing.performance %>% 
  group_by(participant_id) %>% 
  dplyr::summarise(
    n = n(),
    mean_int_accuracy = mean(max_abs_interval_error)
  ) %>% 
  left_join(data_exp1_demographics, by = c("participant_id"="participant_id"))


# exp10
data_exp10_demographics =
  data_exp10 %>% 
  select(participant_id, age, gender, GMSI_MT, GMSI_SA, mt) %>% 
  distinct(participant_id, .keep_all = T)

exp10_batch1 <- read_csv("data/singing-performance-test/experiment10/exp10_batch1.csv")
exp10_batch2 <- read_csv("data/singing-performance-test/experiment10/exp10_batch2.csv")
exp10_batch3 <- read_csv("data/singing-performance-test/experiment10/exp10_batch3.csv")

exp10_batch1_sing = prepare_accuracy_data(exp10_batch1, "b1")
exp10_batch2_sing = prepare_accuracy_data(exp10_batch2, "b2")
exp10_batch3_sing = prepare_accuracy_data(exp10_batch3, "b3")

exp10_sing.performance = as_tibble(rbind(exp10_batch1_sing, exp10_batch2_sing, exp10_batch3_sing))
unique(exp10_sing.performance[, c("participant_id")])  # 201

get_singing_perforamnce_scores(exp10_sing.performance)


exp10_sing.performance_cor =  exp10_sing.performance %>% 
  group_by(participant_id) %>% 
  dplyr::summarise(
    n = n(),
    mean_int_accuracy = mean(max_abs_interval_error)
  ) %>% 
  left_join(data_exp10_demographics, by = c("participant_id"="participant_id"))


# exp11
data_exp11_demographics =
  data_exp11 %>% 
  select(participant_id, age, gender, GMSI_MT, GMSI_SA, mt) %>% 
  distinct(participant_id, .keep_all = T) %>% 
  mutate(participant_id = paste0(participant_id, "b1"))

exp11_batch1 <- read_csv("data/singing-performance-test/experiment11/exp11_batch1.csv")

exp11_sing.performance = prepare_accuracy_data(exp11_batch1, "b1")

unique(exp11_sing.performance[, c("participant_id")])  # 50

get_singing_perforamnce_scores(exp11_sing.performance)


exp11_sing.performance_cor =  exp11_sing.performance %>% 
  group_by(participant_id) %>% 
  dplyr::summarise(
    n = n(),
    mean_int_accuracy = mean(max_abs_interval_error)
  ) %>% 
  left_join(data_exp11_demographics, by = c("participant_id"="participant_id"))

