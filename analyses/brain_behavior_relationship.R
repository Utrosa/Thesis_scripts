library(readxl)
library(openxlsx)
library(ggplot2)
set.seed(123)
## ---------------------------- DEFINE FUNCTIONS ---------------------------- ##
# Compute Spearman's and Pearson's correlation between FC and questionnaire data
compute_correlations <- function(static_fiq, scores, group_name) {
  results <- data.frame(
    Connection = 1:nrow(static_fiq),
    Pearson_r = numeric(nrow(static_fiq)),
    Pearson_p = numeric(nrow(static_fiq)),
    Spearman_rho = numeric(nrow(static_fiq)),
    Spearman_p = numeric(nrow(static_fiq)),
    Group = group_name
  )
  
  for (i in 1:nrow(static_fiq)) {
    static_fc <- as.numeric(static_fiq[i, ])
    pearson_test <- cor.test(static_fc, scores, method = "pearson")
    spearman_test <- cor.test(static_fc, scores, method = "spearman")
    
    results$Pearson_r[i] <- pearson_test$estimate
    results$Pearson_p[i] <- pearson_test$p.value
    results$Spearman_rho[i] <- spearman_test$estimate
    results$Spearman_p[i] <- spearman_test$p.value
  }
  
  return(results)
}

# Check normality of scores
check_normality <- function(data) {
  for (i in 1:nrow(data)) {
    static_fc <- as.numeric(data[i, ])
    # Plot histogram
    hist(static_fc, main = paste("Histogram of Region", i), xlab = "Correlation Coefficient")
    # Q-Q Plot
    qqnorm(static_fc, main = paste("Q-Q Plot of Region", i))
    qqline(static_fc)
    # Shapiro-Wilk Test
    shapiro_result <- shapiro.test(static_fc)
    print(paste("Shapiro-Wilk test for Region", i, ": W =", 
                shapiro_result$statistic,", p-value =",shapiro_result$p.value))
  }
}

# FDR correction for a list of data sets
Padjust <- function(data_frames, p_value_column, method = "BH") {
  # Adjust p-values and add them as new columns
  adjusted_data_frames <- lapply(data_frames, function(df) {
    adjusted_p_values <- p.adjust(df[[p_value_column]], method = method)
    df$adjusted_p_values <- adjusted_p_values
    return(df)
  })
  
  # Combine all data frames into one
  combined_data_frame <- do.call(rbind, adjusted_data_frames)
  
  return(combined_data_frame)
}
# Plot correlation results with FDR correction
plot_spearman_results <- function(data, group_label,title_suffix,color_name,
                                  x_axis_name, x_labels = NULL,
                                  save_path = NULL, save_format = "png") {

  # Adjust p-values for FDR (BH method)
  data$Spearman_p_adjusted <- p.adjust(data$Spearman_p, method = "fdr")
  
  # Create the plot
  plot <- ggplot(data, aes(x = factor(Connection), y = Spearman_rho)) +
    geom_point(aes(color = Spearman_p_adjusted < 0.05), size = 3) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = color_name),
                       breaks = c(FALSE, TRUE),
                       labels = c("Not Significant", "Significant"),
                       name = "FDR Correction") +
    labs(title = paste("Relationship Between FC and", 
                       title_suffix),
         x = x_axis_name,
         y = "Spearman's ρ") +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(size = 20, hjust = 0, margin = margin(b = 20)),
          axis.title.x = element_text(size = 18, margin = margin(t = 15)),
          axis.title.y = element_text(size = 18, margin = margin(t = 15)),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    coord_cartesian(ylim = c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = "dotted")
    
    # Add custom x-labels
    if (!is.null(x_labels)) {
      plot <- plot + scale_x_discrete(labels = x_labels)
    }
  
    # Print or save the plot
    if (!is.null(save_path)) {
      ggsave(
        file.path(save_path, paste0("plot_", group_label, "_", title_suffix, ".", save_format)),
        plot,
        device = save_format,
        width = 8,
        height = 6
      )
    } else {
      print(plot)
    }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GET THE DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
#Get the phenotypic data from the selected participants
datadir = "C:/Users/monik/Documents/Paris/Academic Year (2023-24)/M2 Internship/"
setwd(datadir)
collections <- c("GU_1", "KKI_1", "OHSU_1")
for (collection_name in collections) {
  path <- sprintf("%s/", collection_name)
  path <- paste0(datadir, "Data-collection/", path)
  setwd(path)
  filename <- sprintf("%s_matched_participants.tsv", collection_name)
  assign(collection_name, read.table(filename,
                                     sep='\t',
                                     strip.white=TRUE,
                                     header=TRUE,
                                     fill=TRUE))
}

#Join the pheno data frames into one
all_pheno=rbind(GU_1,KKI_1,OHSU_1)

#Subset the data based on subject type: autistic patients vs healthy controls
#For diagnostic group: 1 = Autism; 2 = Control.
subset_hs=(all_pheno$dx_group==2)
subset_asd=(all_pheno$dx_group==1)

#Create new data frames per subject type
ABIDE_TD=all_pheno[subset_hs,]
ABIDE_ASD=all_pheno[subset_asd,]

## STATIC FC DATA --------------------------------------------------------------
FC_dir = "Analysis/final_results/Static_FC/"
FC_path_ASD = paste0(datadir, FC_dir, "staticFC_subASD_unique.xlsx")
FC_path_TD = paste0(datadir, FC_dir, "staticFC_subTD_unique.xlsx")
static_FC_ASD <- read_xlsx(FC_path_ASD)
static_FC_TD <- read_xlsx(FC_path_TD)

### Check normality of correlation coefficients for all 11 static edges
check_normality(static_FC_ASD[,colnames(static_FC_ASD) %in% ABIDE_ASD$bids_name])
check_normality(static_FC_TD[, colnames(static_FC_TD) %in%  ABIDE_TD$bids_name])

## DYNAMIC FC DATA --------------------------------------------------------------
qpp_dir = "Analysis/final_results/Dynamic_FC/"
qpp_path_ASD = paste0(datadir, qpp_dir, "QPP_metrics_ASD.xlsx")
qpp_path_TD = paste0(datadir, qpp_dir, "QPP_metrics_TD.xlsx")
dynamic_FC_ASD <- read_xlsx(qpp_path_ASD)
dynamic_FC_TD <- read_xlsx(qpp_path_TD)

### Check normality of QPP metrics (strength & frequency)
hist(dynamic_FC_TD$strengthCs)
shapiro.test(dynamic_FC_TD$strengthCs)
qqnorm(dynamic_FC_TD$strengthCs)
qqline(dynamic_FC_TD$strengthCs)

hist(dynamic_FC_TD$frequencyCs)
shapiro.test(dynamic_FC_TD$frequencyCs)
qqnorm(dynamic_FC_TD$frequencyCs)
qqline(dynamic_FC_TD$frequencyCs)

hist(dynamic_FC_ASD$strengthCs)
shapiro.test(dynamic_FC_ASD$strengthCs)
qqnorm(dynamic_FC_ASD$strengthCs)
qqline(dynamic_FC_ASD$strengthCs)

hist(dynamic_FC_ASD$frequencyCs)
shapiro.test(dynamic_FC_ASD$frequencyCs)
qqnorm(dynamic_FC_ASD$frequencyCs)
qqline(dynamic_FC_ASD$frequencyCs)

### Prepare dynamic data frames for analysis
# ASD participants
dynamic_FC_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% 
                                   c("strengthCs", "frequencyCs","subjectID")]
dynamic_FC_ASD <- t(dynamic_FC_ASD)
colnames(dynamic_FC_ASD) <- dynamic_FC_ASD[nrow(dynamic_FC_ASD), ]
dynamic_FC_ASD <- dynamic_FC_ASD[-nrow(dynamic_FC_ASD), ]

# Control participants
dynamic_FC_TD <- dynamic_FC_TD[, colnames(dynamic_FC_TD) %in% 
                                   c("strengthCs", "frequencyCs","subjectID")]
dynamic_FC_TD <- t(dynamic_FC_TD)
colnames(dynamic_FC_TD) <- dynamic_FC_TD[nrow(dynamic_FC_TD), ]
dynamic_FC_TD <- dynamic_FC_TD[-nrow(dynamic_FC_TD), ]

## QUESTIONNAIRE SCORES --------------------------------------------------------
### INTELIGENCE ----------------------------------------------------------------
#### Check which questionnaires were used to measure full-scale IQ (FIQ)
cat(unique(GU_1$fiq_test_type))
cat(unique(KKI_1$fiq_test_type))
cat(unique(OHSU_1$fiq_test_type))

#### Count ASD and TD participants for the common questionnaire: WISC-IV
#### The minimal sample size per group is 21. dx_group: 1=Autism; 2=Control 
gu_fiq1 <- sum(GU_1$fiq_test_type == "WISC-IV" & GU_1$dx_group == 1)
kki_fiq1 <- sum(KKI_1$fiq_test_type == "WISC-IV" & KKI_1$dx_group == 1)
ohsu_fiq1 <- sum(OHSU_1$fiq_test_type == "WISC-IV" & OHSU_1$dx_group == 1)

gu_fiq2 <- sum(GU_1$fiq_test_type == "WISC-IV" & GU_1$dx_group == 2)
kki_fiq2 <-sum(KKI_1$fiq_test_type == "WISC-IV" & KKI_1$dx_group == 2)
ohsu_fiq2 <- sum(OHSU_1$fiq_test_type == "WISC-IV" & OHSU_1$dx_group == 2)

cat(sum(gu_fiq1,kki_fiq1,ohsu_fiq1), sum(gu_fiq2,kki_fiq2,ohsu_fiq2))

#### Extract the participant IDs and scores for valid (comparable) FIQ data
asd_fiq <- ABIDE_ASD$bids_name[ABIDE_ASD$fiq_test_type == "WISC-IV"]
td_fiq <-  ABIDE_TD$bids_name[ABIDE_TD$fiq_test_type == "WISC-IV"]
asd_fiq_vec <- ABIDE_ASD$fiq[ABIDE_ASD$fiq_test_type == "WISC-IV"]
td_fiq_vec <- ABIDE_TD$fiq[ABIDE_TD$fiq_test_type == "WISC-IV"]

#### Subset static_FC for the selected participants (comparable questionnaires)
static_fiq_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% asd_fiq]
static_fiq_TD <- static_FC_TD[, colnames(static_FC_TD) %in% td_fiq]

#### Subset dynamic_FC for the selected participants (comparable questionnaires)
dynamic_fiq_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% asd_fiq]
dynamic_fiq_TD <- dynamic_FC_TD[, colnames(dynamic_FC_TD) %in% td_fiq]

#### Check for normality of fiq data: non-sig. == normality.
shapiro.test(asd_fiq_vec)
shapiro.test(td_fiq_vec)
qqnorm(asd_fiq_vec)
qqline(asd_fiq_vec)
qqnorm(td_fiq_vec)
qqline(td_fiq_vec)

#### Compute correlation between static FC and FIQ scores for all 11 connections
results_fiq_TD <- compute_correlations(static_fiq_TD, td_fiq_vec, "TD")
results_fiq_ASD <- compute_correlations(static_fiq_ASD, asd_fiq_vec, "ASD")
combined_results_fiq <- rbind(results_fiq_TD, results_fiq_ASD)

#### Compute correlation between dynamic FC and FIQ scores for all 11 connections
results_fiq_TD2 <- compute_correlations(dynamic_fiq_TD, td_fiq_vec, "TD")
results_fiq_ASD2 <- compute_correlations(dynamic_fiq_ASD, asd_fiq_vec, "ASD")
combined_results_fiq2 <- rbind(results_fiq_TD2, results_fiq_ASD2)

#### Save the results
write.xlsx(combined_results_fiq, file = paste0(datadir,FC_dir,
                                               "static_fc_fiq_correlations.xlsx"))
write.xlsx(combined_results_fiq2, file = paste0(datadir,qpp_dir,
                                               "dynamic_fc_fiq_correlations.xlsx"))
### ADI-R ----------------------------------------------------------------------
#### Extract ADI scores (ASD participants only)
adi_a_asd <- as.numeric(ABIDE_ASD$adi_r_social_total_a)
adi_verbal_asd <- as.numeric(ABIDE_ASD$adi_r_verbal_total_bv)
adi_nonverbal_asd <- as.numeric(ABIDE_ASD$adi_r_nonverbal_total_bv)
adi_c_asd <- as.numeric(ABIDE_ASD$adi_r_rrb_total_c)
adi_d_asd <- as.numeric(ABIDE_ASD$adi_r_onset_total_d)

#### Check for normality of ADI scores
shapiro.test(adi_a_asd)
shapiro.test(adi_verbal_asd)
shapiro.test(adi_nonverbal_asd)
shapiro.test(adi_c_asd)
shapiro.test(adi_d_asd)

#### Subset static & dynamic FC for all ASD participants (the same questionnaire)
asd_all <- ABIDE_ASD$bids_name
static_adi_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% asd_all]
dynamic_adi_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% asd_all]

#### Compute correlation between static FC and ADI scores for all 11 connections
results_ADIa <- compute_correlations(static_adi_ASD, adi_a_asd, "ASD")
results_ADIbv <- compute_correlations(static_adi_ASD, adi_verbal_asd, "ASD")
results_ADIbnv <- compute_correlations(static_adi_ASD, adi_nonverbal_asd, "ASD")
results_ADIc <- compute_correlations(static_adi_ASD, adi_c_asd, "ASD")
results_ADId <- compute_correlations(static_adi_ASD, adi_d_asd, "ASD")

#### Apply FDR correction
combined_results_ADI <- Padjust(list(results_ADIa, results_ADIbv, results_ADIbnv,
             results_ADIc, results_ADId),"Spearman_p",method="BH")
head(combined_results_ADI)

#### Compute correlation between dynamic FC and ADI scores for all 11 connections
results_ADIa2 <- compute_correlations(dynamic_adi_ASD, adi_a_asd, "ASD")
results_ADIbv2 <- compute_correlations(dynamic_adi_ASD, adi_verbal_asd, "ASD")
results_ADIbnv2 <- compute_correlations(dynamic_adi_ASD, adi_nonverbal_asd, "ASD")
results_ADIc2 <- compute_correlations(dynamic_adi_ASD, adi_c_asd, "ASD")
results_ADId2 <- compute_correlations(dynamic_adi_ASD, adi_d_asd, "ASD")
combined_results_ADI2 <- rbind(results_ADIa2, results_ADIbv2, results_ADIbnv2,
                               results_ADIc2, results_ADId2)

#### Apply FDR correction & save the adjusted p-value to the table
combined_results_ADI2 <- Padjust(list(results_ADIa2, results_ADIbv2, results_ADIbnv2,
             results_ADIc2, results_ADId2),"Spearman_p",method="BH")
head(combined_results_ADI2)

#### Save the results
write.xlsx(combined_results_ADI,file=paste0(datadir,FC_dir,
                                            "static_fc_ADI_correlations.xlsx"))
write.xlsx(combined_results_ADI2,file=paste0(datadir,qpp_dir,
                                            "dynamic_fc_ADI_correlations.xlsx"))
### ADOS -----------------------------------------------------------------------
#### Was ADOS scored/administered by research reliable personnel? 0=no; 1=yes.
cat(unique(GU_1$ados_rsrch_reliable[GU_1$dx_group == 1]))
cat(unique(KKI_1$ados_rsrch_reliable[GU_1$dx_group == 1]))
cat(unique(OHSU_1$ados_rsrch_reliable[GU_1$dx_group == 1]))

#### Which ADOS module was used to assess ASD symptoms?
cat(unique(GU_1$ados_module))
cat(unique(KKI_1$ados_module))
cat(unique(OHSU_1$ados_module))

#### Extract the participant IDs and scores for which there's ADOS-2 data
ABIDE_ASD$ados_2_total[ABIDE_ASD$ados_2_total == "n/a"] <- 9999
ABIDE_ASD$ados_2_total <- as.numeric(ABIDE_ASD$ados_2_total)
ados2 <- ABIDE_ASD$bids_name[ABIDE_ASD$ados_2_total != 9999 & 
                             ABIDE_ASD$ados_module == 3]
ados2_scores <- ABIDE_ASD$ados_2_total[ABIDE_ASD$ados_2_total != 9999 &
                                       ABIDE_ASD$ados_module == 3]

#### Extract the participant IDs and scores for which there's ADOS-G data
ABIDE_ASD$ados_g_total[ABIDE_ASD$ados_g_total == "n/a"] <- 9999
ABIDE_ASD$ados_g_total <- as.numeric(ABIDE_ASD$ados_g_total)
adosG <- ABIDE_ASD$bids_name[ABIDE_ASD$ados_g_total != 9999 &
                             ABIDE_ASD$ados_module == 3]
adosG_scores <- ABIDE_ASD$ados_g_total[ABIDE_ASD$ados_g_total != 9999 &
                                       ABIDE_ASD$ados_module == 3]
#### Ascertain that n > 21
length(ados2) > 21
length(adosG) > 21

#### Subset static & dynamic FC for valid scores (the same ADOS for all)
static_ados2_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% ados2]
static_adosG_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% adosG]
dynamic_ados2_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% ados2]
dynamic_adosG_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% adosG]

#### Check for normality of ADOS scores
shapiro.test(ados2_scores)
shapiro.test(adosG_scores)
qqnorm(ados2_scores)
qqline(ados2_scores)
qqnorm(adosG_scores)
qqline(adosG_scores)

#### Compute correlation between static FC and ADOS scores for all 11 edges
results_ados2 <- compute_correlations(static_ados2_ASD, ados2_scores, "ASD")
results_adosG <- compute_correlations(static_adosG_ASD, adosG_scores, "ASD")
combined_results_ados <- Padjust(list(results_ados2, results_adosG),"Spearman_p",method="BH")

#### Compute correlation between QPP and ADOS scores for all 11 edges
results_ados22 <- compute_correlations(dynamic_ados2_ASD, ados2_scores, "ASD")
results_adosG2 <- compute_correlations(dynamic_adosG_ASD, adosG_scores, "ASD")
combined_results_ados2 <- Padjust(list(results_ados22, results_adosG2),"Spearman_p",method="BH")

#### Save the results
write.xlsx(combined_results_ados, file = paste0(datadir,FC_dir,
                                                "static_fc_ados_correlations.xlsx"))
write.xlsx(combined_results_ados2, file = paste0(datadir,qpp_dir,
                                               "dynamic_fc_ados_correlations.xlsx"))

### SRS ------------------------------------------------------------------------
#### Which SRS version was used? 1=child from; 2=adult form
cat(unique(GU_1$srs_version))
cat(unique(KKI_1$srs_version))
cat(unique(OHSU_1$srs_version))

#### Which SRS edition was used?
cat(unique(GU_1$srs_edition))
cat(unique(KKI_1$srs_edition))
cat(unique(OHSU_1$srs_edition))

#### Extract the participant IDs and scores for which there's SRS-1 data
srs1 <- ABIDE_ASD$bids_name[ABIDE_ASD$srs_edition == 1 &
                                ABIDE_ASD$srs_version == 1]
srs2 <- ABIDE_ASD$bids_name[ABIDE_ASD$srs_edition == 2 &
                                ABIDE_ASD$srs_version == 1]

#### Extract the participant IDs and scores for which there's SRS-2 data
srs1_scores <- ABIDE_ASD$srs_total_raw[ABIDE_ASD$srs_edition == 1 &
                                       ABIDE_ASD$srs_version == 1]
srs2_scores <- ABIDE_ASD$srs_total_raw[ABIDE_ASD$srs_edition == 2 &
                                       ABIDE_ASD$srs_version == 1]
srs1_scores <- as.numeric(srs1_scores)
srs2_scores <- as.numeric(srs2_scores)

#### Ascertain that n > 21
length(srs1) > 21
length(srs2) > 21

#### Subset static FC & QPP for valid SRS scores (the same questionnaire)
static_srs1_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% srs1]
static_srs2_ASD <- static_FC_ASD[, colnames(static_FC_ASD) %in% srs2]
dynamic_srs1_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% srs1]
dynamic_srs2_ASD <- dynamic_FC_ASD[, colnames(dynamic_FC_ASD) %in% srs2]

#### Check for normality of SRS scores
shapiro.test(srs1_scores)
shapiro.test(srs2_scores)
qqnorm(srs1_scores)
qqline(srs1_scores)
qqnorm(srs2_scores)
qqline(srs2_scores)

#### Compute correlation between static FC and SRS scores for all 11 edges
results_srs1 <- compute_correlations(static_srs1_ASD, srs1_scores, "ASD")
results_srs2 <- compute_correlations(static_srs2_ASD, srs2_scores, "ASD")
combined_results_srs <- Padjust(list(results_srs1, results_srs2),"Spearman_p",method="BH")

#### Compute correlation between dynamic FC and SRS scores for all 11 edges
results_srs12 <- compute_correlations(dynamic_srs1_ASD, srs1_scores, "ASD")
results_srs22 <- compute_correlations(dynamic_srs2_ASD, srs2_scores, "ASD")
combined_results_srs2 <- Padjust(list(results_srs12, results_srs22),"Spearman_p",method="BH")

#### Save the results
write.xlsx(combined_results_srs, file = paste0(datadir,FC_dir,
                                              "static_fc_srs_correlations.xlsx"))
write.xlsx(combined_results_srs2, file = paste0(datadir,qpp_dir,
                                               "dynamic_fc_srs_correlations.xlsx"))
## PLOTTING --------------------------------------------------------------------
# To add x-axis labels to static FC plots:
## static_ASD_FC$pairs_net7
## theme(axis.text.x = element_text(angle = 45, hjust = 1)
FC_dir = "Analysis/final_results/Static_FC/brain-behavior"
qpp_dir = "Analysis/final_results/Dynamic_FC/brain-behavior/"

plot_spearman_results(results_fiq_TD,"TD","FIQ (TD)","red","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_fiq_ASD,"ASD","FIQ (ASD)","blue","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))

plot_spearman_results(results_fiq_TD2,"TD","FIQ (TD)","red","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_fiq_ASD2,"ASD","FIQ (ASD)","blue","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))

plot_spearman_results(results_srs1,"SRS1","SRS1","green","Static Functional Connection Index", x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_srs2,"SRS2","SRS2", "green","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))

plot_spearman_results(results_srs12,"SRS1","SRS1","green","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_srs22,"SRS2","SRS-2", "green","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))

plot_spearman_results(results_ados2,"ADOS2","ADOS2","orange","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_adosG,"ADOSG","ADOSG","orange","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))

plot_spearman_results(results_ados22,"ADOS2","ADOS2","orange","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_adosG2,"ADOSG","ADOSG","orange","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))

plot_spearman_results(results_ADIa,"ADIa","ADIa","purple","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_ADIbv,"ADIbv","ADIb verbal","purple","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_ADIbnv,"ADIbnv","ADIb non-verbal","purple","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_ADIc,"ADIc","ADIc","purple","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))
plot_spearman_results(results_ADId,"ADId","ADId","purple","Static Functional Connection Index",x_labels = NULL,paste0(datadir,FC_dir))

plot_spearman_results(results_ADIa2,"ADIa","ADIa","purple","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_ADIbv2,"ADIbv","ADIb verbal","purple","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_ADIbnv2,"ADIbnv","ADI-R non-verbal","purple","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_ADIc2,"ADIc","ADIc","purple","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))
plot_spearman_results(results_ADId2,"ADId","ADId","purple","QPP Metrics",c("QPP Strength","QPP Frequency"),paste0(datadir,qpp_dir))