################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Supporting methods for melodic features extraction via bootstrap
################################################################################
sum_boot.data_to.plot = function(data, vars_to_sum){
  data_sum = data %>%
    group_by(degree) %>%
    dplyr::summarise_at(vars(vars_to_sum), list(mean = mean, sd = sd), na.rm = TRUE)
  return(data_sum)
}


get_bootstrapped_features <- function(
    data_singing, 
    data_melcon,
    vars_sung_intervals,
    nboot, 
    # bw,
    # interval_range,  # only activate if using the kde in peak finding
    is_free_notes = FALSE
){
  
  degrees = table(data_singing$degree)
  
  if (is_free_notes) {
    data.ready = data_singing %>%
      select(participant_id:degree, network_id, list_sung_ints, 
             root_mean_squared_interval)
  } else {
    data.ready = data_singing %>%
      unite("list_sung_ints", vars_sung_intervals, sep=",", remove=FALSE) %>% 
      select(participant_id:degree, network_id, list_sung_ints, 
             vars_sung_intervals, root_mean_squared_interval)
  }

  store <- c()
  
  # loop over iterations
  for (i in 1:nboot){
    
    print(paste("features boot", i, "out of", nboot))
    
    # sample with replacement
    data_sample <- data.ready %>% group_by(network_id) %>% group_split()
    data_sample <- sample(data_sample, length(data_sample), replace=TRUE)
    data_sample <- bind_rows(data_sample) 
    
    if (is_free_notes) {
      data_sample.long = data_sample %>%
        mutate(row_id = row_number()) %>% 
        unnest(list_sung_ints) %>% 
        rename(interval = list_sung_ints)
    } else {
      data_sample.long = data_sample %>%
        mutate(row_id = row_number()) %>% 
        pivot_longer(vars_sung_intervals, names_to = "index", values_to = "interval") 
    }

    # loop over generations
    
    # we calculate peaks in matlab now
    # boot.peaks_by_degree = c()
    
    # entropy
    boot.entropies_by_degree = c()
    
    degrees = table(data.ready$degree)
    
    for (d in 1:length(degrees)){
      
      degree.d = d-1
      print(paste0("degree ", degree.d, ", in boot ", i))
      
      data_sample.long_degree = data_sample.long %>%  filter(degree == degree.d)

      # peaks
      
      # smooth
      # data_sample.long_degree_kde <- density(data_sample.long_degree[["interval"]],
      #                                        bw = bw,
      #                                        kernel = "gaussian",
      #                                        from = min(interval_range) - 3*bw,
      #                                        to = max(interval_range) + 3*bw
      # )
      
      
      # boot.peaks_by_degree[[d]] = get_peaks_iterations(data_sample.long_degree_kde)
      
      # this method computes entropies at different granularities (but is slow)
      # boot.entropies_by_degree[[d]] = get_entropy_iterations(data_sample.long_degree)
      
      boot.entropies_by_degree[[d]] = get_entropy_0.25(data_sample.long_degree)
      
    }
    
    # peaks_table = do.call(rbind, boot.peaks_by_degree)
    entropies = do.call(rbind, boot.entropies_by_degree)
    
    # melodic consonance
    data_sample.long_round = data_sample.long  %>%
      mutate(interval = round(interval,1)) %>% 
      select(degree, interval)
    
    data_sum_melcon = as_tibble(
      merge(data_sample.long_round, data_melcon, by = "interval")) %>%
      group_by(degree) %>%dplyr::summarise(m = mean(mean_rating))
    
    data_int_size = data_sample.long %>% 
      group_by(row_id) %>% 
      mutate(mean_abs_int_size = mean(abs(interval))) %>% 
      group_by(degree) %>%
      dplyr::summarise(
        abs_int_size = mean(mean_abs_int_size, na.rm =T)
      )
    
    data_int_error = data_sample %>%
      group_by(degree) %>%
      dplyr::summarise(
        RMSE_interval = mean(root_mean_squared_interval, na.rm =T)
      )
    
    store[[i]] = tibble(
      boot = i,
      degree = data_int_error$degree,
      # RMSE
      RMSE_interval = data_int_error$RMSE_interval,
      # interval size
      abs_int_size = data_int_size$abs_int_size,
      # melodic pleasantness
      melcon = data_sum_melcon$m,
      # peaks
      # n_peaks=peaks_table$n_peaks,
      # boot.peaks=peaks_table$peaks,
      # values=peaks_table$values,
      # entropy
      entropy_0.25=entropies,
    )
    
  }
  o = do.call(rbind, store)
  return(o)
}


get_bootstrapped_features_simulations <- function(
    data, 
    vars_sung_intervals,
    nboot
){
  
  degrees = table(data$degree)
  
  store <- c()
  
  # loop over iterations
  for (i in 1:nboot){
    
    print(paste("features boot", i, "out of", nboot))
    
    # sample with replacement
    data_sample <- data %>% group_by(melody_id) %>% group_split()
    data_sample <- sample(data_sample, length(data_sample), replace=TRUE)
    data_sample <- bind_rows(data_sample) 
    
    data_sample.long = data_sample %>%
      mutate(row_id = row_number()) %>% 
      pivot_longer(vars_sung_intervals, names_to = "index", values_to = "interval") 
    
    # entropy
    boot.entropies_by_degree = c()
    
    degrees = table(data$degree)
    
    for (d in 1:length(degrees)){
      
      degree.d = d-1
      print(paste0("degree ", degree.d, ", in boot ", i))
      
      data_sample.long_degree = data_sample.long %>%  filter(degree == degree.d)
      
      boot.entropies_by_degree[[d]] = get_entropy_0.25(data_sample.long_degree)
      
    }
    
    entropies = do.call(rbind, boot.entropies_by_degree)
    
    # size
    data_int_size = data_sample.long %>% 
      group_by(row_id) %>% 
      mutate(mean_abs_int_size = mean(abs(interval))) %>% 
      group_by(degree) %>%
      dplyr::summarise(
        abs_int_size = mean(mean_abs_int_size, na.rm =T)
      )
    
    # error
    data_int_error = data_sample %>%
      group_by(degree) %>%
      dplyr::summarise(
        RMSE_interval = mean(rmse_interval, na.rm =T)
      ) %>% 
      mutate(RMSE_interval = ifelse(degree==0,NA,RMSE_interval))
    
    store[[i]] = tibble(
      boot = i,
      degree = data_int_error$degree,
      # RMSE
      RMSE_interval = data_int_error$RMSE_interval,
      # interval size
      abs_int_size = data_int_size$abs_int_size,
      # entropy
      entropy_0.25=entropies
    )
    
  }
  o = do.call(rbind, store)
  return(o)
}


get_proportion_large_ints_boot <- function(
    data_singing, 
    vars_sung_intervals,
    nboot
){
  
  degrees = table(data_singing$degree)
  
  data.ready = data_singing %>%
    unite("list_sung_ints", vars_sung_intervals, sep=",", remove=FALSE) %>% 
    select(participant_id:degree, network_id, list_sung_ints, 
           vars_sung_intervals, root_mean_squared_interval)
  
  store <- c()
  
  # loop over iterations
  for (i in 1:nboot){
    
    print(paste("features boot", i, "out of", nboot))
    
    # sample with replacement
    data_sample <- data.ready %>% group_by(network_id) %>% group_split()
    data_sample <- sample(data_sample, length(data_sample), replace=TRUE)
    data_sample <- bind_rows(data_sample) 
    
    store[[i]]  = data_sample %>% 
      pivot_longer(vars_sung_intervals)  %>% 
      mutate(abs_int = abs(value)) %>% 
      mutate(is_larger.than.7 = ifelse(abs_int > 7, 1 , 0)) %>% 
      group_by(degree) %>% 
      dplyr::summarise(
        n_large.7 = sum(is_larger.than.7), 
        n = n(),
        proportion = (n_large.7 / n) * 100
      )
    
  }
  o = do.call(rbind, store)
  return(o)
}


get_peaks_iterations <- function(data){
  
  peaks = get_num_peaks(data$x, data$y)
  
  x_points = c()
  y_points = c()
  for (j in 1:length(peaks[[1]])){
    x = peaks[[1]][[j]]$interval
    y = peaks[[1]][[j]]$value
    
    x_points = c(x_points, x)
    y_points = c(y_points, y)
    
    peaks_df = tibble(x = x_points, y = y_points) 
    
    n_peaks = length(peaks_df$x)
    
    peaks_table = tibble(
      i=j, 
      n_peaks=n_peaks, 
      peaks=list(peaks_df$x), 
      values=list(peaks_df$y)
    )
    
  }
  
  return(peaks_table)
}


get_entropy_iterations = function(data){
  
  intervals = data %>% 
    select(degree, interval)
  
  
  H.int.0 = c()
  H.int.0.25 = c()
  H.int.0.5 = c()
  H.int.0.1 = c()
  H.continous = c()
  
  # entropy
  dat_h_int_0 = round(intervals$interval)
  dat_h_int_0.25 = round(intervals$interval/0.25)*0.25
  dat_h_int_0.5 = round(intervals$interval/0.5)*0.5
  dat_h_int_0.1 = round(intervals$interval, 1)
  
  c.h_int_0 = categorical_entropy(dat_h_int_0)
  c.h_int_0.25 = categorical_entropy(dat_h_int_0.25)
  c.h_int_0.5 = categorical_entropy(dat_h_int_0.5)
  c.h_int_0.1 = categorical_entropy(dat_h_int_0.1)
  
  H.int.0 = c(H.int.0, c.h_int_0)
  H.int.0.25 = c(H.int.0.25, c.h_int_0.25)
  H.int.0.5 = c(H.int.0.5, c.h_int_0.5)
  H.int.0.1 = c(H.int.0.1, c.h_int_0.1)
  
  # continous method
  cont.entropies = entropy(intervals$interval, k = 10)
  cont.entropies_clean = cont.entropies[!is.na(cont.entropies) & !is.infinite(cont.entropies)]
  mean_cont.entropies = mean(cont.entropies_clean, na.rm = T)
  
  H.continous = c(H.continous, mean_cont.entropies)
  
  # store entropies
  entropies = tibble(
    H.int.0 = H.int.0,
    H.int.0.5 = H.int.0.5,
    H.int.0.25 = H.int.0.25,
    H.int.0.1 = H.int.0.1,
    H.continous = H.continous
  )
  
  return(entropies)
}


get_entropy_0.25 = function(data){
  
  intervals = data %>% 
    select(degree, interval)
  
  
  H.int.0.25 = c()

  # entropy
  dat_h_int_0.25 = round(intervals$interval/0.25)*0.25

  c.h_int_0.25 = categorical_entropy(dat_h_int_0.25)

  H.int.0.25 = c(H.int.0.25, c.h_int_0.25)
  
  return(H.int.0.25)
}


categorical_entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}


get_features_bootstraped_slider <- function(
    data_singing,
    data_melcon,
    vars_sung_intervals,
    nboot
    # bw,
    # interval_range
){
  
  degrees = table(data_singing$degree)
  
  data.ready = data_singing %>%
    select(participant_id:degree, network_id, 
           vars_sung_intervals, interval_error)
  
  store <- c()
  
  # loop over iterations
  for (i in 1:nboot){
    
    print(paste("features boot", i, "out of", nboot))
    
    # sample with replacement
    data_sample <- data.ready %>% group_by(network_id) %>% group_split()
    data_sample <- sample(data_sample, length(data_sample), replace=TRUE)
    data_sample <- bind_rows(data_sample) 
    
    # loop over generations
    # boot.peaks_by_degree = c()
    boot.entropies_by_degree = c()
    
    for (d in 1:length(degrees)){
      
      degree.d = d-1
      print(paste0("degree ", degree.d, ", in boot ", i))
      
      data_sample_degree = data_sample %>%  filter(degree == degree.d)

      # smooth
      # data_sample.degree_kde <- density(data_sample_degree[[vars_sung_intervals]],
      #                                        bw = bw,
      #                                        kernel = "gaussian",
      #                                        from = min(interval_range) - 3*bw,
      #                                        to = max(interval_range) + 3*bw
      # )
      # 
      # # peaks
      # boot.peaks_by_degree[[d]] = get_peaks_iterations(data_sample.degree_kde)
      
      # entropies
      boot.entropies_by_degree[[d]] = get_entropy_0.25(data_sample_degree)
      
    }
    
    # peaks_table = do.call(rbind, boot.peaks_by_degree)
    entropies = do.call(rbind, boot.entropies_by_degree)
    
    
  # melcon
    # melodic consonance
    data_sample_round = data_sample  %>%
      mutate(interval = round(interval,1)) %>% 
      select(degree, interval)
    
    data_sum_melcon = as_tibble(
      merge(data_sample_round, data_melcon, by = "interval")) %>%
      group_by(degree) %>%dplyr::summarise(m = mean(mean_rating))
    
    # others
    scores = data_sample %>%  
      group_by(degree) %>%
      dplyr::summarise(
        mean_abs_int_size = mean(abs(interval), na.rm =T),
        mean_error_interval = mean(interval_error, na.rm =T)
      )
    
    store[[i]] = tibble(
      boot = i,
      degree = scores$degree,
      # RMSE
      mean_error_interval = scores$mean_error_interval,
      # interval size
      mean_abs_int_size = scores$mean_abs_int_size,
      # peaks
      # peaks=peaks_table$n_peaks,
      # melcon
      melcon = data_sum_melcon$m,
      # information theory
      interval_entropy_0.25=entropies[,1]
    )
    
  }
  o = do.call(rbind, store)
  return(o)
}


any_input_list_to_numeric = function(input_list){
  if (is.character(input_list)){
    numeric_list = convert_string_numbers_to_numeric_vector(input_list)
  } else if (is.list(input_list)) {
    numeric_list = unlist(input_list)
  } else {
    numeric_list = input_list
  }
  return(numeric_list)
}

calculate_mean_abs_int_size = function(intervals_raw){
  ints = any_input_list_to_numeric(intervals_raw)
  actual_intervals = abs(ints)
  mean_abs_intervals = mean(actual_intervals)
  return(mean_abs_intervals)
}

se <- function(x) sqrt(var(x) / length(x))

convert_string_numbers_to_numeric_vector = function(intervals_string){
  library(comprehenr)
  intervals_string_sep = as.vector(strsplit(intervals_string, ",")[[1]]) # convert to strings separated
  list_intervals_num = to_vec(for(i in intervals_string_sep) as.numeric(i)) # convert to numeric
  return(list_intervals_num)
}


compute_error <- function(raw_diffs) {
  # Calculate the absolute differences
  abs_differences <- abs(raw_diffs)
  
  # Calculate the maximum error
  max_error <- ifelse(length(abs_differences) == 0, 999, max(abs_differences))
  
  # Calculate the root mean squared error
  mean_squared_differences <- ifelse(length(abs_differences) == 0, 999, mean(abs_differences^2))
  root_mean_squared <- ifelse(length(abs_differences) == 0, 999, sqrt(mean_squared_differences))
  
  # Return the maximum error and root mean squared error
  return(root_mean_squared)
}


