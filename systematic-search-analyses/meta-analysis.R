######### ------------------------------------------------------------ #########
#                                                                              #
#             R script for funnel and forest plot for effect sizes             #
#                                                                              #
######### ------------------------------------------------------------ #########

# Load the library for meta-analysis
library(metafor)
library(dplyr)

#Read in the data file. This file needs to be in tab-separated value form.
df <- read.csv("Meta-Analysis-Static-FC.tsv",sep="\t")

# See an overview of the data
summary(df)

# Create two dataframes: static FC & static FC - phenotype correlation
df_static <- subset(df,effect_type == 'FC_static')
df_pheno <- subset(df,effect_type == 'FC-phenotype')

#### --------------------- Estimate Outcome Measures ---------------------- ####
# The outcome measures are yi (effect size) and vi (sampling variance).
# Note that escalc() corrects yi for positive bias, yielding Hedges' g.

# Contrast two independent groups
res_static <- escalc(measure = "SMD", 
              di = as.numeric(cohen_d), 
              n1i = participant_group_N,
              n2i = control_N, 
              data = df_static,
              vtype='AV') # AV OR UB aren't defaults

# Describe the direction and strength of the association between two variables
res_pheno <- escalc(measure = "COR", # Raw correlation coeff.
                    ri = as.numeric(raw_corr),
                    ti = as.numeric(ti_pheno),
                    ni = sample_size,
                    data = df_pheno,
                    vtype='AV') # AV OR UB aren't defaults

#### ------------------- Model Fitting: Random Effects -------------------- ####
RE_static = rma(yi = res_static$yi,
                vi = res_static$vi,
                data = res_static,
                slab = res_static$paper_shortname)

RE_pheno = rma(yi = res_pheno$yi,
               vi = res_pheno$vi,
               data = res_pheno,
               slab = res_pheno$paper_shortname)

#### ---------------------- Forest and Funnel Plots ----------------------- ####
png("forest_RE_static.png", width = 2500, height = 3500, res = 300) 
op <- par(cex = 0.75, font = 2)
forest(RE_static,
       header="Author(s) and Year",
       shade=TRUE)
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(op)
dev.off()

png("forest_RE_pheno.png", width = 6500, height = 8500, res = 300) 
op <- par(cex = 1, font = 2)
forest(RE_pheno,
       header="Author(s) and Year",
       shade=TRUE)
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(op)
dev.off()

png("funnel_RE_static.png", width = 1500, height = 1500, res = 300) 
funnel(RE_static)
dev.off()

png("funnel_RE_pheno.png", width = 1500, height = 1500, res = 300)
funnel(RE_pheno)
dev.off()


# Tests for funnel plot asymmetry
regtest(RE_static)
regtest(RE_pheno)

#### ------------ Forest & Funnel plots: Mean Effect per Study ------------ ####

# Aggregate multiple effect sizes within the same study ------------------------
res_studies_static <- aggregate(x = res_static, 
                               cluster = paper_shortname,
                               struct="CS",
                               rho=0.6)

res_studies_pheno <- aggregate(x = res_pheno, 
                              cluster = paper_shortname,
                              struct="CS",
                              rho=0.6)
# Note. Since outcomes within clusters cannot be assumed to be independent, the 
# struct parameter is set to "CS". An arbitrary 0.6 value is used for rho.

# Model fitting for aggregated effect sizes ------------------------------------
RE_static_mean = rma(yi = res_studies_static$yi,
                     vi = res_studies_static$vi,
                     data = res_studies_static,
                     slab = res_studies_static$paper_shortname)

RE_pheno_mean = rma(yi = res_studies_pheno$yi,
                    vi = res_studies_pheno$vi,
                    data = res_studies_pheno,
                    slab = res_studies_pheno$paper_shortname)

# Plotting for aggregated effect sizes -----------------------------------------
png("forest_RE_static_mean.png", width = 1500, height = 1500, res = 300) 
op <- par(cex = 0.8, font = 2)
forest(RE_static_mean,
       header="Author(s) and Year",
       shade=TRUE)
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(op)
dev.off()

png("forest_RE_pheno_mean.png", width = 2500, height = 2500, res = 300) 
op <- par(cex = 1,5, font = 2)
forest(RE_pheno_mean,
       header="Author(s) and Year",
       shade=TRUE)
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(op)
dev.off()

png("funnel_RE_static_mean.png", width = 1500, height = 1500, res = 300) 
funnel(RE_static_mean)
dev.off()

png("funnel_RE_pheno_mean.png", width = 1500, height = 1500, res = 300)
funnel(RE_pheno_mean)
dev.off()
###############################################################################



