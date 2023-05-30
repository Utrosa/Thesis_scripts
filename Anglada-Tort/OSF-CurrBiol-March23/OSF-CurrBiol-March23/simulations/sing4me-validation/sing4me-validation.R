################################################################################
# Current Biology (2022)
# Authors: Manuel Anglada-Tort, Peter Harrison, and Nori Jacoby
# Script: Analysis sing4me validation
################################################################################


################################################################################
# Global
################################################################################
# imports
library(tidyverse)
library(egg)
library(ggpubr)

set.seed(2022)
loadNamespace("egg")
theme_set(theme_pubr())

# methods
current_path = "simulations/sing4me-validation/"


################################################################################
# plot output validation
################################################################################
sing4me_val <- read_csv("simulations/sing4me-validation/sing4me_simulations.csv")

mean(abs(sing4me_val$accuracy))
sd(sing4me_val$accuracy)
min(sing4me_val$accuracy)
max(sing4me_val$accuracy)

plot = ggplot(sing4me_val, aes(x=target_tone, y=accuracy)) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(40, 70, 2)) +
  ylim(-0.025, 0.025) +
  geom_vline(xintercept = seq(40,70,5), colour = "lightgrey", linetype="dashed", size = 0.25) +
  geom_hline(yintercept = 0, colour = "darkred", linetype="dashed", size = 0.5) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())  

ggsave(paste0(current_path, "simulations_sing4m.pdf"),  
       plot, 
       height = 8, width = 15, units = "cm", dpi=300)
