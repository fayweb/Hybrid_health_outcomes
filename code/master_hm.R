# ***********************************************************
# Title: Predicting the health outcomes of infections in hybrid mice

# Purpose: This script defines all the settings and executes
#         all the code (.R, .md) to reproduce the analysis
#         of the project
#
# Authors: Fay Webster
# ID variables:
# ***********************************************************
# Part 1: Set standard settings & install packages            ----
# ***********************************************************
  # Install packages/load libraries to keep R environment stable
    # install
        # pacman for simplified bulk pkg import
        # renv for pkg consistency over time
#install.packages("pacman")
# increase maximum overlaps
options(ggrepel.max.overlaps = Inf)

library(pacman)
    ## Standard settings ----
        # seed
set.seed(13102023)

        # Use p_load to install (if not already) and load the packages
pacman::p_load(mice, stringr, gridExtra, dplyr, tidyverse, tidyr, janitor, 
               visdat, corrplot, RColorBrewer, ggplot2, VIM, limma, 
               latticeExtra, patchwork,FactoMineR, ggrepel, factoextra, 
               reshape2, sjPlot, stargazer, jtools,modelsummary, ggeffects, 
               pheatmap, ggpubr, ggridges, gt, caret, randomForest, rfUtilities)
    
## Define within project file paths ----
        # code
c <- "code"
clab      <- paste0(c, "/lab/")
cfield     <- paste0(c, "/field/")
canalysis <- paste0(c, "/analysis/")
cdesign <- paste0(c, "/design/") # experimental project design
nmi   <- paste0(c, "/nmi/")
cmodels <- paste0(c, "/models/")

        # data
            # building dynamic paths
                # Get the user's profile directory on Windows
user_profile <- Sys.getenv("USERPROFILE")

                # Append the specific path
one_drive <- file.path(user_profile, "OneDrive",
                       "Documents", "GitHub", "Hybrid_health_outcomes")

                # relative_path is the desired path
d <- paste0(one_drive, "/data")

            # labs
dlab <- paste0(d, "/lab")
dlab_raw <- paste0(dlab, "/raw")
  dlab_inter <- paste0(dlab, "/intermediate")
  dlab_final <- paste0(dlab, "/final")

            # field
dfield <- paste0(d, "/field")
  dfield_raw <- paste0(dfield, "/raw")
  dfield_inter <- paste0(dfield, "/intermediate")
  dfield_final <- paste0(dfield, "/final")


            # data product for analysis
danalysis <- paste0(d, "/analysis")
    danal_final <- paste0(danalysis, "/final") 

        # output
output <- paste0(one_drive, "/output")
  figures <- paste0(output, "/figures")
    fi <- paste0(figures, "/imputation")
    an_fi <- paste0(figures, "/analysis")
    d_fi <- paste0(figures, "/design")
    panels_fi <- paste0(figures, "/panels")
  tables  <- paste0(output, "/tables")
  
  
 #### vectors for selecting genes for analysis
 Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
                 "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                 "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                 "TICAM1", "TNF")

    ## Define functions ----
if (1) source(file.path(c, "functions.R"))

# ***********************************************************
# Part 2: Run Data cleaning - lab                        ----
# ***********************************************************
#----------------------------------------------------------*
# 2.1: Import raw data & save as intermediate/processed
#----------------------------------------------------------*
if (1) source(file.path(clab, "lab_import.R"))
#----------------------------------------------------------*
# 2.2: Conduct cleaning (formatting) 
# Creates: Challenge
#----------------------------------------------------------*
if (1) source(file.path(clab, "lab_clean.R"))
#----------------------------------------------------------*
# 2.3: Visualize data
 #----------------------------------------------------------*
if (1) source(file.path(clab, "lab_visualize.R"))
# Creates: Correlation matrix between laboratory gene expression values
# ***********************************************************
# Part 3: Run field infection data cleaning                      ----
# **********************************************************
#----------------------------------------------------------*
# 3.1: Import raw data & save as intermediate/processed
# Requires:
# Creates: field_cleaned_data.csv
#----------------------------------------------------------*
if (1) source(file.path(cfield, "field_import.R"))
#----------------------------------------------------------*
# 3.2: Conduct cleaning (formatting) w/o changing data
#----------------------------------------------------------*
if (1) source(file.path(cfield, "field_clean.R"))
 #----------------------------------------------------------*
# 3.3: Visualize data
#----------------------------------------------------------*
 if (1) source(file.path(cfield, "field_visualize.R"))
 #----------------------------------------------------------*
 # 3.4: Import amplicon infection intensities and join with field
 #----------------------------------------------------------*
 if (1) source(file.path(cfield, "field_import_amplicon_intensities.R"))

 
# ***********************************************************
# Part 4:  MNI: Merge normalize impute                           ----
# ***********************************************************
#----------------------------------------------------------*
# 4.1 Merging of field and laboratory data
 # Creates output: merge_prior_imputation.csv
#----------------------------------------------------------*
if (1) source(file.path(nmi, "nmi_merge_long.R"))
# 4.2 Create a dataframe with the selection of the mice in the experiments
# Creates output: mice_selection.csv
#----------------------------------------------------------*
if (1) source(file.path(nmi, "nmi_merge_wide.R"))
# 4.3 Normalization of data 
# Creates output: genes
#----------------------------------------------------------*
if (1) source(file.path(nmi, "nmi_normalize.R"))
# 4.4 Imputation of missing data 
# Creates output: genes
#----------------------------------------------------------*
if (1) source(file.path(nmi, "nmi_impute.R"))
# 4.4 Imputation of missing data 
# Creates output: hm
#----------------------------------------------------------*
 
 
 ######### Decision: Removing IL.10 from gene selection due data largely missing
 # remove il_10
 #### vectors for selecting genes for analysis
 Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF")
 
 # ***********************************************************
 # Part 5: Experimental design                           ----
 # ***********************************************************
 #----------------------------------------------------------*
 # Show the primary results of our experimental design
 # How many rodents, distributions, strains, and parasite information
 if (0) source(file.path(cdesign, "design_experimental.R"))
 
 
# ***********************************************************
# Part 6: Analysis                           ----
# ***********************************************************
#----------------------------------------------------------*
# 6.1: PCA 
 # performing a pca analysis on the laboratory immune gene data
# Requires: hm
# Creates: lab
 # Plots: biplot, pca_variables,
 #  contr_PC1, contr_PC2
#----------------------------------------------------------*
if (0) source(file.path(canalysis, "analysis_PCA_genes_lab.R"))
 # 6.2: PCA 
 # Regressions with pc axes 
 # Plots: pc1_current_infection, pc2_current_infection, coefs5
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, "analysis_linear_regressions_PCA.R"))
#----------------------------------------------------------*
# 6.2: Heatmap lab genes
# Requires: hm, lab
#----------------------------------------------------------*
if (0) source(file.path(canalysis, "heatmap_lab_genes.R"))
 #----------------------------------------------------------*
 # 6.3: Multiple multivariate regression of genes vs weight loss in the lab
 # Requires: hm, lab
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, "analysis_multiple_multivariate_regression.R"))
 

 
 # ***********************************************************
 # Part 7: Analysis                           ----
 # ***********************************************************
 #----------------------------------------------------------*
 #7: PCA 
 # performing a pca analysis on the laboratory immune gene data
 # Requires: hm
 # Creates: lab
 # Plots: biplot, pca_variables,
 #  contr_PC1, contr_PC2
 #----------------------------------------------------------*
 

# ***********************************************************
# Part end: Save information key information about project    ----
# ***********************************************************
sessionInfo()
