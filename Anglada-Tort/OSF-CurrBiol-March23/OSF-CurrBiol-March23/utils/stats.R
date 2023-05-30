################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: General methods for stats
################################################################################


################################################################################
# stats features bootstrap
################################################################################
run_lm_bootraped_data  = function(data, column){
  
  store = c()
  
  nboot =  length(unique(data$boot))
  
  for (i in 1:nboot){
    
    data_sample = data %>% filter(boot == i)
    
    lm.entropy = lm(entropy_0.25 ~as.numeric(degree), data=data_sample)
    b.entropy = lm.entropy$coefficients[[2]]
    
    lm.n_peaks = lm(n_peaks ~as.numeric(degree), data=data_sample)
    b.n_peaks = lm.n_peaks$coefficients[[2]]
    
    lm.abs_int_size = lm(abs_int_size ~as.numeric(degree), data=data_sample)
    b.abs_int_size = lm.abs_int_size$coefficients[[2]]
    
    lm.RMSE_interval = lm(RMSE_interval ~as.numeric(degree), data=data_sample)
    b.RMSE_interval = lm.RMSE_interval$coefficients[[2]]
    
    store[[i]] = tibble(
      boot = i,
      b.entropy = b.entropy,
      b.n_peaks = b.n_peaks,
      b.abs_int_size = b.abs_int_size,
      b.RMSE_interval = b.RMSE_interval
    )
  }
  bs = do.call(rbind, store)
  
  
  stats = as_tibble(
    rbind(CI_boot(bs$b.entropy),
          CI_boot(bs$b.n_peaks),
          CI_boot(bs$b.abs_int_size),
          CI_boot(bs$b.RMSE_interval))
  ) %>% 
    mutate(feture = c("entropy", "peaks", "int.size", "error"))
  
  return(list(bs, stats))
}


CI_boot = function(x){
  tibble(
    mean = mean(x),
    se  =  sd(x),
    lower_95.CI  = mean - se * 1.96,
    upper_95.CI  = mean + se * 1.96
  )
}


################################################################################
# functions to power analysis
################################################################################
run_sample.size_analysis = function(data, bw, N_list, N_boot){
  
  data_power_marginals = data %>% 
    filter(degree %in% 8:10) %>% 
    select(participant_id, network_id, degree, sung_interval1, sung_interval2) 
  
  data_power_features = data %>% 
    select(participant_id, network_id, degree, sung_interval1, sung_interval2)
  
  store = c()
  store2 = c()
  
  for (i in 1:N_boot){
    for (n in 1:length(N_list)){
      
      N =  N_list[n]
      
      print(paste0("iteration ", i , " with N = ", N))
      
      # 1. slice data based on N  
      data_power_marginals_sliced <- data_power_marginals %>% 
        nest(cols = -network_id) %>% 
        slice_sample(n=N) %>% 
        unnest(cols = -network_id) %>% 
        pivot_longer(sung_interval1:sung_interval2, names_to = "int", values_to = "interval") %>% 
        select(-int)
      
      data_power_features_sliced <- data_power_features %>% 
        nest(cols = -network_id) %>% 
        slice_sample(n=N) %>% 
        unnest(cols = -network_id) %>% 
        pivot_longer(sung_interval1:sung_interval2, names_to = "int", values_to = "interval") %>% 
        select(-int)
      
      # 2. split chains by half
      
      # marginals
      sample_size_marginals = floor(0.5*nrow(table(data_power_marginals_sliced$network_id)))
      rand_selected_marginals = sample(names(table(data_power_marginals_sliced$network_id)), size = sample_size_marginals)
      
      half1_marginals <- data_power_marginals_sliced %>% filter(network_id %in% rand_selected_marginals) %>%  drop_na()
      half2_marginals <- data_power_marginals_sliced %>% filter(!network_id %in% rand_selected_marginals) %>%  drop_na()
      
      # features
      sample_size_features = floor(0.5*nrow(table(data_power_features_sliced$network_id)))
      rand_selected_features = sample(names(table(data_power_features_sliced$network_id)), size = sample_size_features)
      
      half1_features <- data_power_features_sliced %>% filter(network_id %in% rand_selected_features) %>%  drop_na()
      half2_features <- data_power_features_sliced %>% filter(!network_id %in% rand_selected_features) %>%  drop_na()
      
      # 3. Entropy correlation
      degrees = table(half1_features$degree)
      entropies_by_degree_half1 = c()
      entropies_by_degree_half2 = c()
      
      for (d in 1:length(degrees)){
        
        degree.d = d-1
        print(paste0("degree ", degree.d, ", in boot ", i))
        
        half1_features_degree = half1_features %>%  filter(degree == degree.d)
        half2_features_degree = half2_features %>%  filter(degree == degree.d)
        
        # entropy
        entropies_by_degree_half1[[d]] = get_entropy_0.25(half1_features_degree)
        entropies_by_degree_half2[[d]] = get_entropy_0.25(half2_features_degree)
      }
      
      total_entropies_by_degree_half1 = do.call(rbind, entropies_by_degree_half1)
      total_entropies_by_degree_half2 = do.call(rbind, entropies_by_degree_half2)
      
      entropy.cor.r = cor.test(total_entropies_by_degree_half1[,1], total_entropies_by_degree_half2[,1])
      entropy.cor.sp = cor.test(total_entropies_by_degree_half1[,1], total_entropies_by_degree_half2[,1], method = "spearman")
      
      # 4. Smooth  data
      print(paste("smoothing iteration", i, ", N=", N, "..."))
      print(paste("half1 N chains = ",  length(table(half1_marginals$network_id))))
      print(paste("half2 N chains = ",  length(table(half2_marginals$network_id))))
      
      half1.kde <- density(half1_marginals[["interval"]],
                           bw = bw,
                           kernel = "gaussian",
                           from = min(interval_range) - 3*bw,
                           to = max(interval_range) + 3*bw)
      
      half1.smoothed <- data.frame(
        x = half1.kde$x,
        y = half1.kde$y
      )
      
      half2.kde <- density(half2_marginals[["interval"]],
                           bw = bw,
                           kernel = "gaussian",
                           from = min(interval_range) - 3*bw,
                           to = max(interval_range) + 3*bw)
      
      half2.smoothed <- data.frame(
        x = half2.kde$x,
        y = half2.kde$y
      )
      
      # 5. Peaks only correlation
      print(paste("peakscorrelation, iteration", i, ", N=", N, "..."))
      
      gam.half1 <- gam(y ~ s(x), data = half1.smoothed, method = "REML")
      gam.half2 <- gam(y ~ s(x), data = half2.smoothed, method = "REML")
      
      half1.smoothed$curve = gam.half1$fitted.values
      half1.smoothed$rest = half1.smoothed$y - half1.smoothed$curve
      
      half2.smoothed$curve = gam.half2$fitted.values
      half2.smoothed$rest = half2.smoothed$y - half2.smoothed$curve
      
      cor.res.gam = tibble(
        var1=half1.smoothed$rest,
        var2=half2.smoothed$rest
      )
      
      split.cor.res.gam = item_split_half(cor.res.gam)
      
      # 6. JSD
      print(paste("jsd, iteration", i, ", N=", N, "..."))
      jsd = run_jsd(half1.smoothed, half2.smoothed)
      
      store[[n]] = tibble(
        N= N,
        boot=i,
        entropy.cor.r = entropy.cor.r$estimate,
        entropy.cor.sp = entropy.cor.sp$estimate,
        peaks_splithalf = split.cor.res.gam$splithalf,
        peaks_spearmanbrown = split.cor.res.gam$spearmanbrown,
        jsd = jsd[[1]]
      )
      
      output1 = do.call(rbind, store)
    }
    store2[[i]] = tibble(
      output = output1,
      i = i
    )
  }
  
  res = do.call(rbind, store2)
  
  res$N = res$output$N
  res$boot = res$output$boot
  res$entropy.cor.r = res$output$entropy.cor.r
  res$entropy.cor.sp = res$output$entropy.cor.sp
  res$peaks_splithalf = res$output$peaks_splithalf
  res$peaks_spearmanbrown = res$output$peaks_spearmanbrown
  res$jsd = res$output$jsd
  
  return(res)
}


summarize_results_boot = function(data){
  results_power_sum = data %>% 
    group_by(N) %>% 
    dplyr::summarise(
      n = n(),
      # mean
      JSD=mean(output$jsd),
      entropy.cor.r=mean(output$entropy.cor.r),
      entropy.cor.sp=mean(output$entropy.cor.sp),
      peaks_splithalf=mean(output$peaks_splithalf),
      peaks_spearmanbrown=mean(output$peaks_spearmanbrown)
    ) %>% 
    select(N, JSD:peaks_spearmanbrown) %>% 
    pivot_longer(JSD:peaks_spearmanbrown, "metric", values_to = "mean") 
  
  results_power_sum$metric = factor(
    results_power_sum$metric, 
    levels = c("JSD", "peaks_splithalf", "peaks_spearmanbrown", "entropy.cor.r", "entropy.cor.sp"))
  
  return(results_power_sum)
}


# jsd
run_jsd = function(d1, d2){
  P <- d1$y/sum(d1$y)
  Q <- d2$y/sum(d2$y)
  PQ <- rbind(P,Q)
  jsd <- JSD(PQ)
  return(jsd)
}


plot_power.analysis_metric = function(data, metric_str, color_bar, y_lab, title){
  plot_jsd = data %>% 
    filter(metric == metric_str) %>% 
    ggplot(aes(x=N, y=mean, group = metric, fill=metric)) + 
    geom_bar(stat="identity", color="black", fill=color_bar, position=position_dodge()) +
    ylab(y_lab) +
    xlab("N chains") +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(size=8))
  
  return(plot_jsd)
}


# Melodic pleasantness ratings
smooth_melcon_ratings = function(data, n_bootstrap=1000){
  dyads_by_participant <- data %>% 
    group_by(participant_id) %>% group_split()
  
  dyads_smoothed_bootstrapped <- map_dfr(1:n_bootstrap, function(i) {
    n_participants <- length(dyads_by_participant)
    ind <- sample(n_participants, n_participants, replace = TRUE)
    data_bootstrapped <- bind_rows(dyads_by_participant[ind])
    tibble(
      iteration = i,
      interval =seq(min(interval_range),max(interval_range),by=0.1),
      rating = smooth_1d_gaussian(data_bootstrapped$target_interval,data_bootstrapped$answer,interval,BW)
    )
  })
  
  dyads_smoothed_final <- 
    dyads_smoothed_bootstrapped %>% 
    group_by(interval) %>% 
    dplyr::summarise(
      mean_rating = mean(rating, na.rm = TRUE),
      se_rating = sd(rating, na.rm = TRUE)
    )
  return(dyads_smoothed_final)
}


run_sample.size_analysis_ratings = function(data, bw, N_list, N_boot){
  
  store = c()
  store2 = c()
  
  for (i in 1:N_boot){
    for (n in 1:length(N_list)){
      
      N =  N_list[n]
      
      print(paste0("iteration ", i , " with N = ", N))
      
      # 1. slice data based on N  
      data_power_sliced <- data_exp6 %>% 
        slice_sample(n=N) 
      
      # 2. split chains by half
      
      # marginals
      sample_size_marginals = floor(0.5*nrow(table(data_power_sliced$stimulus_id)))
      rand_selected_marginals = sample(names(table(data_power_sliced$stimulus_id)), size = sample_size_marginals)
      
      half1_marginals <- data_power_sliced %>% filter(stimulus_id %in% rand_selected_marginals) %>%  drop_na()
      half2_marginals <- data_power_sliced %>% filter(!stimulus_id %in% rand_selected_marginals) %>%  drop_na()
      
      # 3. Smooth  data
      print(paste("smoothing iteration", i, ", N=", N, "..."))
      print(paste("half1 N chains = ",  length(table(half1_marginals$stimulus_id))))
      print(paste("half2 N chains = ",  length(table(half2_marginals$stimulus_id))))
      
      y = seq(min(interval_range),max(interval_range),by=0.1)
      
      half1.kde = smooth_1d_gaussian(
        half1_marginals$target_interval,
        half1_marginals$answer,
        y, bw)
      
      half2.kde = smooth_1d_gaussian(
        half2_marginals$target_interval,
        half2_marginals$answer,
        y, bw)
      
      half1.smoothed <- data.frame(
        x = y,
        y = half1.kde
      )
      
      half2.smoothed <- data.frame(
        x = y,
        y = half2.kde
      )
      
      # 5. Peaks only correlation
      print(paste("peakscorrelation, iteration", i, ", N=", N, "..."))
      
      gam.half1 <- gam(y ~ s(x), data = half1.smoothed, method = "REML")
      gam.half2 <- gam(y ~ s(x), data = half2.smoothed, method = "REML")
      
      half1.smoothed$curve = gam.half1$fitted.values
      half1.smoothed$rest = half1.smoothed$y - half1.smoothed$curve
      
      half2.smoothed$curve = gam.half2$fitted.values
      half2.smoothed$rest = half2.smoothed$y - half2.smoothed$curve
      
      cor.res.gam = tibble(
        var1=half1.smoothed$rest,
        var2=half2.smoothed$rest
      )
      
      split.cor.res.gam = item_split_half(cor.res.gam)
      
      # 6. JSD
      print(paste("jsd, iteration", i, ", N=", N, "..."))
      jsd = run_jsd(half1.smoothed, half2.smoothed)
      
      store[[n]] = tibble(
        N= N,
        boot=i,
        peaks_splithalf = split.cor.res.gam$splithalf,
        peaks_spearmanbrown = split.cor.res.gam$spearmanbrown,
        jsd = jsd[[1]]
      )
      
      output1 = do.call(rbind, store)
    }
    store2[[i]] = tibble(
      output = output1,
      i = i
    )
  }
  
  res = do.call(rbind, store2)
  
  res$N = res$output$N
  res$boot = res$output$boot
  res$peaks_splithalf = res$output$peaks_splithalf
  res$peaks_spearmanbrown = res$output$peaks_spearmanbrown
  res$jsd = res$output$jsd
  
  return(res)
}


################################################################################
# functions to power analysis
################################################################################
run_jsd_comparison = function(data1, data2, bw, interval_range, nboot){
  
  data1_ready = data1 %>% 
    filter(degree  %in% 8:10) %>% 
    select(network_id, sung_interval1:sung_interval2) %>%  
    pivot_longer(sung_interval1:sung_interval2, values_to = "interval")
  
  data2_ready = data2 %>% 
    filter(degree  %in% 8:10) %>% 
    select(network_id, sung_interval1:sung_interval2) %>%  
    pivot_longer(sung_interval1:sung_interval2, values_to = "interval")
  
  store = c()
  
  # boot
  for (i in 1:nboot){
    
    print(paste("boot", i, "out of", nboot))
    
    # sample with replacement
    data1_sample <- data1_ready %>% group_by(network_id) %>% group_split()
    data1_sample <- sample(data1_sample, length(data1_sample), replace=TRUE)
    data1_sample <- bind_rows(data1_sample) 
    
    data2_sample <- data2_ready %>% group_by(network_id) %>% group_split()
    data2_sample <- sample(data2_sample, length(data2_sample), replace=TRUE)
    data2_sample <- bind_rows(data2_sample) 
    
    
    # smooth
    print(paste("smoothing iteration", i))
    print(paste("half1 N chains = ",  length(table(data1_sample$network_id))))
    print(paste("half2 N chains = ",  length(table(data2_sample$network_id))))
    
    sample1.kde <- density(data1_sample[["interval"]],
                         bw = bw,
                         kernel = "gaussian",
                         from = min(interval_range) - 3*bw,
                         to = max(interval_range) + 3*bw)
    
    sample1.smoothed <- data.frame(
      x = sample1.kde$x,
      y = sample1.kde$y
    )
    
    sample2.kde <- density(data2_sample[["interval"]],
                         bw = bw,
                         kernel = "gaussian",
                         from = min(interval_range) - 3*bw,
                         to = max(interval_range) + 3*bw)
    
    sample2.smoothed <- data.frame(
      x = sample2.kde$x,
      y = sample2.kde$y
    )
    
    # JSD
    print(paste("jsd, iteration", i))
    jsd = run_jsd(sample1.smoothed, sample2.smoothed)
    
    store[[i]] = jsd[[1]]
  }
  output = do.call(rbind, store)
  return(output)
}
