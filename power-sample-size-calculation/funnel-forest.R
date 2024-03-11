######### ------------------------------------------------------------ #########
#                                                                              #
#             R script for funnel and forest plot for effect sizes             #
#                                                                              #
######### ------------------------------------------------------------ #########

# Load the library for meta-analysis
library(metafor)

#Read in the data file. This file needs to be in tab-separated value form.
df <- read.csv("Meta-Analysis-Static-FC.tsv",sep="\t")

# See an overview of the data
summary(df)

#### ------------------ Forest & Funnel plots: All Data ------------------- ####

# Extract the required variables: yi (effect size), vi (standard error).
yi_static <- as.numeric(df$cohen_d[df$effect_type == "FC_static"])
vi_static <- 1 / sqrt(df$sample_size[df$effect_type == "FC_static"])

yi_pheno <- as.numeric(df$cohen_d[df$effect_type == "FC-phenotype"])
vi_pheno <- 1 / sqrt(df$sample_size[df$effect_type == "FC-phenotype"])

# Static Functional Connectivity
forest(yi_static, vi_static,
       slab=(df$paper_shortname[df$effect_type == "FC_static"]))
title(main = "Static FC (all)")

funnel(yi_static, vi_static)
title(main = "Static FC (all)")

# Correlation Between Static Functional Connectivity and Clinical Symptoms
forest(yi_pheno, vi_pheno,
       slab=df$paper_shortname[df$effect_type == "FC-phenotype"])
title(main = "Correlation Between Static FC and Clinical Symptoms (all)")

funnel(yi_pheno, vi_pheno)
title(main = "Correlation Between Static FC and Clinical Symptoms (all)")

#### ------------- Forest & Funnel plots: Mean Data per Study ------------- ####

# Extract the required variables: yi (effect size), vi (standard error).
yi_static <- as.numeric(df$cohen_d[df$effect_type == "FC_static"])
vi_static <- 1 / sqrt(df$sample_size[df$effect_type == "FC_static"])

yi_pheno <- as.numeric(df$cohen_d[df$effect_type == "FC-phenotype"])
vi_pheno <- 1 / sqrt(df$sample_size[df$effect_type == "FC-phenotype"])

# Calculate mean values per category of $paper_shortname
mean_yi_static <- tapply(yi_static, 
                         df$paper_shortname[df$effect_type == "FC_static"],
                         mean)
mean_vi_static <- tapply(vi_static, 
                         df$paper_shortname[df$effect_type == "FC_static"],
                         mean)

mean_yi_pheno <- tapply(yi_pheno, 
                         df$paper_shortname[df$effect_type == "FC-phenotype"],
                         mean)
mean_vi_pheno <- tapply(vi_pheno,
                        df$paper_shortname[df$effect_type == "FC-phenotype"],
                        mean)
#### ----------- Create Forest & Funnel plots ----------- ####

# Static Functional Connectivity
forest(mean_yi_static, mean_vi_static,
       slab=(unique(df$paper_shortname[df$effect_type == "FC_static"])),
       header="Author(s) and Year", shade=TRUE)
title(main = "Static FC (mean)")

funnel(mean_yi_static, mean_vi_static)
title(main = "Static FC (mean)")

# Correlation Between Static Functional Connectivity and Clinical Symptoms
forest(mean_yi_pheno, mean_vi_pheno,
       slab=unique(df$paper_shortname[df$effect_type == "FC-phenotype"]),
       header="Author(s) and Year", shade = TRUE)
title(main = "Correlation Between Static FC and Clinical Symptoms (mean)")

funnel(mean_yi_pheno, mean_vi_pheno)
title(main = "Correlation Between Static FC and Clinical Symptoms (mean)")

###############################################################################



