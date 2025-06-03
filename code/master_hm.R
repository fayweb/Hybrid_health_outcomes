# ***********************************************************
# Title: Predicting the Health Outcomes of Parasite Infections in Hybrid Mice
#
# Purpose: This master script initializes the project environment,
#          including all standard settings, package installations,
#          data paths, and custom functions required to conduct
#          comprehensive analyses. It sets a consistent and reproducible
#          foundation for importing, cleaning, visualizing, normalizing,
#          imputing, analyzing, and modeling infection health outcomes
#          in laboratory and wild hybrid mouse populations.
#
# Workflow Structure:
#   1. Standard settings: Set seeds, load libraries, define paths
#   2. Custom function definitions: Create functions to aid visualization
#      and statistical distribution testing.
#   3. Prepare data paths: Dynamically define file paths for efficient
#      and reproducible data handling.
#   4. Data cleaning and preparation (lab and field)
#   5. Merge, normalize, and impute data
#   6. Laboratory infection analysis (PCA, linear models)
#   7. Random forest model development and validation
#   8. Wild mouse analysis and cross-population validation
#
# Author: Fay Webster
# Date: Initiated October 13, 2023
# ***********************************************************

# ***********************************************************
# Part 1: Set Standard Settings & Load Packages ----
# ***********************************************************

# Install packages/load libraries to maintain a stable R environment
options(ggrepel.max.overlaps = Inf)

library(pacman)

# Seed for reproducibility
set.seed(13102023)

# Load necessary packages using pacman
pacman::p_load(mice, stringr, gridExtra, dplyr, tidyverse, tidyr, janitor, 
               visdat, corrplot, RColorBrewer, ggplot2, VIM, limma, 
               latticeExtra, patchwork, FactoMineR, ggrepel, factoextra, 
               reshape2, sjPlot, stargazer, jtools, modelsummary, ggeffects, 
               pheatmap, ggpubr, ggridges, gt, caret, randomForest, rfUtilities,
               parasiteLoad, fitdistrplus, optimx, leaflet, magick, ggdist,
               ggbeeswarm, ggtext, kableExtra, webshot, broom, flextable,
               viridis)

# ***********************************************************
# Part 2: Define Project File Paths ----
# ***********************************************************

# Code directories
c <- "code"
clab <- paste0(c, "/lab/")
cfield <- paste0(c, "/field/")
canalysis <- paste0(c, "/analysis/")
cdesign <- paste0(c, "/design/") # Experimental project design
nmi <- paste0(c, "/nmi/")
cmodels <- paste0(c, "/models/")

# Data directories
user_profile <- Sys.getenv("USERPROFILE")
one_drive <- file.path(user_profile, "OneDrive", "Documents", "GitHub", "Hybrid_health_outcomes")
d <- paste0(one_drive, "/data")

# Lab data paths
dlab <- paste0(d, "/lab")
dlab_raw <- paste0(dlab, "/raw")
dlab_inter <- paste0(dlab, "/intermediate")
dlab_final <- paste0(dlab, "/final")

# Field data paths
dfield <- paste0(d, "/field")
dfield_raw <- paste0(dfield, "/raw")
dfield_inter <- paste0(dfield, "/intermediate")
dfield_final <- paste0(dfield, "/final")

# Data product for analysis
danalysis <- paste0(d, "/analysis")
danal_final <- paste0(danalysis, "/final")

# Output directories
output <- paste0(one_drive, "/output")
figures <- paste0(output, "/figures")
fi <- paste0(figures, "/imputation")
an_fi <- paste0(figures, "/analysis")
d_fi <- paste0(figures, "/design")
panels_fi <- paste0(figures, "/panels")
tables <- paste0(output, "/tables")

# Vectors for selecting genes for analysis
Genes_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
             "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
             "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
             "TICAM1", "TNF")

EqPCR.cols <- c("delta_ct_cewe_MminusE", "MC.Eimeria", "Ct.Eimeria")

# ***********************************************************
# Part 3: Define Custom Functions ----
# ***********************************************************
if (1) source(file.path(c, "functions.R"))
# ***********************************************************

# Part 2: Run Data Cleaning - Laboratory Data               ----
# ***********************************************************

#----------------------------------------------------------*
# 2.1: Import raw laboratory data 

#----------------------------------------------------------*
# 2.2: Clean and format lab data
# Creates: Challenge
#----------------------------------------------------------*
if (1) source(file.path(clab, "lab_clean.R"))  # Harmonizes parasite naming, infection histories, etc.

#----------------------------------------------------------*
# 2.3: Visualize gene correlations
#----------------------------------------------------------*
if (1) source(file.path(clab, "lab_visualize.R"))  # Generates gene-gene correlation matrix (mLN only)

# ***********************************************************
# Part 3: Field infection data cleaning and integration      ----
# ***********************************************************
# Purpose: Clean field metadata and enrich it with qPCR 
#          infection intensities and amplicon-based species
#          identifications. Save intermediate and final versions
#          to ensure modular execution.
#
# Requires: Field_infection_data.csv, CEWE_FECES_infection_intensities.txt,
#           Sample_selection_Metabarcoding_Complete.csv
# Creates: field_imported_raw.csv (intermediate)
#          field_cleaned_intermediate.csv (intermediate)
#          field_cleaned_data.csv (final)
#          cor_genes_field.jpeg (figure)
# ***********************************************************

#----------------------------------------------------------*
# 3.1: Import raw data & save as intermediate
#----------------------------------------------------------*
if (1) source(file.path(cfield, "field_import.R"))

#----------------------------------------------------------*
# 3.2: Conduct basic cleaning and formatting
#----------------------------------------------------------*
if (1) source(file.path(cfield, "field_clean.R"))

#----------------------------------------------------------*
# 3.3: Visualize field gene expression correlations
#----------------------------------------------------------*
if (1) source(file.path(cfield, "field_visualize.R"))

#----------------------------------------------------------*
# 3.4: Integrate infection intensities and amplicon species calls
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
 if (1) source(file.path(nmi, "Clean_after_impute.R"))
 # 4.5 Cleaning after imputing  
 # Correct the parasite labels for legend
 # Add custom colours for parasites throughout the scripts
 # Creates factor levels for parasite strains
 #----------------------------------------------------------*

# ***********************************************************
# Part 5: Experimental Design of Laboratory Infections ----
# ***********************************************************
# Purpose: Summarize the experimental setup in the lab infections.
# Outputs: Clean design tables, plots of experimental groups, timepoints

# 5.1: Summarize experimental groups and setup
if (0) source(file.path(cdesign, "design_tables_paper.R"))
if (0) source(file.path(cdesign, "design_experimental.R"))
# 5.2 Measurement Methods: Weight, OPG, qPCR
if (0) source(file.path(cdesign, "design_field_quantification_methods.R"))
# This script explains how weight, oocysts, and parasite burden were quantified
# It can generate a supplementary table or explanatory paragraph

# 5.3 Immune Gene Panel: Target selection and references
if (0) source(file.path(cdesign, "design_immune_genes.R"))
# Details which genes were selected and why
# Can generate a table of target genes, reference genes, primer source, etc.

# 5.4 Eimeria Quantification Pipeline
if (0) source(file.path(cdesign, "design_quantification_eimeria.R"))
# Explains how qPCR or microscopy counts were processed
# Outputs: summary stats, method explanation (can go in methods/supp)

# 5.5 Modeling Strategy: Outcome, predictors, covariates
if (0) source(file.path(cdesign, "design_models.R"))
# Describes the rationale for using weight loss as outcome, infection as predictor
# Sets up linear model structures for downstream analysis (Part 6)


# ***********************************************************
# Part 6: Laboratory Infection Analysis ----
# ***********************************************************
# Purpose: Analyze immune gene expression and health outcomes in lab infections
# This corresponds to Results section "Immune gene expression profiles correlate with infection-induced health impacts"

# 6.1: Individual Gene Responses to Infection
if (0) source(file.path("code/analysis/lab_infections/gene_responses_individual.R"))
# Regression analysis for each gene vs infection status
# Outputs: Table 1, effect sizes, significance levels
# Corresponds to paragraph: "Species-specific immune gene responses to Eimeria infection"

# 6.2: Principal Component Analysis  
if (0) source(file.path("code/analysis/lab_infections/pca_immune_genes.R"))
# PCA on immune gene expression, variance explained, gene loadings
# Outputs: Figure 1B, PC scores, gene contribution tables
# Corresponds to paragraph: "Principal component analysis reveals coordinated immune response patterns"

# 6.3: Linear Models - Immune Signatures Predict Weight Loss
if (0) source(file.path("code/analysis/lab_infections/linear_models_pc_weightloss.R"))
# PC1, PC2 as predictors of weight loss
# Multiple nested models, interaction effects
# Corresponds to paragraph: "Immune signatures predict weight loss in infected mice"
#----------------------------------------------------------*

#----------------------------------------------------------*
# 6.4: Random Forest Model Development
# Train and validate random forest model on laboratory data
# Purpose: Build predictive model for weight loss based on immune gene expression
# Requires: Challenge dataset with complete immune gene data and weight loss outcomes
# Creates: Trained random forest model, performance metrics, variable importance
#----------------------------------------------------------*
if (1) source(file.path("code/analysis/lab_infections/random_forest_training.R"))
#----------------------------------------------------------*
# 6.4.1: Random Forest Model Diagnostics
# Generate diagnostic plots and performance validation for RF model
# Purpose: Validate model assumptions and create supplementary figures
# Requires: Trained RF model, test set predictions, cross-validation results
# Creates: Diagnostic plots (residuals, Q-Q, CV performance), correlation matrix, model comparison table
#----------------------------------------------------------*
if (0) source(file.path("code/analysis/lab_infections/random_forest_diagnostics.R"))
#----------------------------------------------------------*
# 6.5 Random Forest Model Validation in Laboratory Data
# Validate RF predictions against known infection parameters in lab data
# Purpose: Test if RF predictions correlate with infection status, species, and intensity
# Requires: Trained RF model, Challenge dataset
# Creates: Validation plots, correlation statistics, species comparison
#----------------------------------------------------------*
if (0) source(file.path(clab_inf, "random_forest_validation.R"))


# ***********************************************************
# Part 7: Random Forest Application to Wild Mice ----
# ***********************************************************
# Purpose: Apply laboratory-trained RF model to wild mice and validate predictions
# This corresponds to Results section "Random forest predictions in wild-caught mice"
#----------------------------------------------------------*
# 7.1: Apply Random Forest Model to Wild Mice
# Apply lab-trained RF model to predict weight loss in wild mice
# Purpose: Generate weight loss predictions for all wild mice based on immune gene expression
# Requires: Trained RF model (from 6.4), hm dataset (wild mice with immune gene data)
# Creates: hm dataset with predicted weight loss values, application summary stats
#----------------------------------------------------------*
if (1) source(file.path("code/analysis/wild_mice/random_forest_apply_field.R"))
#----------------------------------------------------------*
# 7.2: Validate Predictions with 
# - Infection Intensity
# - Infection status
#----------------------------------------------------------*
if (0) source(file.path("code/analysis/wild_mice/validate_infection_intensity.R"))
#----------------------------------------------------------*
# 7.5: Test Parasite Community Effects
# Analyze whether other gastrointestinal parasites affect RF predictions
# Purpose: Determine if RF predictions are Eimeria-specific or influenced by other parasites
# Requires: hm dataset with RF predictions, complete parasite community data
# Creates: Multi-parasite models, Eimeria-specificity analysis
# Corresponds to: "Eimeria-specific effects in the context of parasite communities"
#----------------------------------------------------------*
if (1) source(file.path("code/analysis/field_application/validate_parasite_specificity.R"))

#----------------------------------------------------------*
# 7.6: Single Gene Validation - CXCL9 Cross-Population Test
# Test whether CXCL9 alone can predict infection costs across populations
# Purpose: Validate the most important RF predictor as a standalone biomarker
# Requires: Challenge dataset (lab), hm dataset (field), CXCL9 expression data
# Creates: Cross-population CXCL9 model, validation statistics, biomarker assessment
# Corresponds to: "CXCL9 as a conserved predictor across populations"
#----------------------------------------------------------*
if (1) source(file.path("code/analysis/field_application/validate_cxcl9_biomarker.R"))


#############
#######                 CAN PROBABLY REMOVE THE NEXT SCRIPTS
###                    Keeping them here until we are close to submission.. 



















 # Now let's analyse the gene expression value distribution
 if (0) source(file = "code/analysis/Analysis_PCA_lab_infections.R")
 # Can the PC1 and PC2 representing immune gene expression predict weight loss?
 if (0) source(file = "code/analysis/Analysis_linear_regressions_PCA_laboratory_infections.R")
 # Let's create the panel figure 1 for the manuscript
 if (0) source(file = "code/figure_creation/Panel_1.R")
 
 
 # ***********************************************************
 # Part 6: How different are wild-derived and wild-caugth mice?
 # ***********************************************************
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, 
                         "analysis_multiple_multivariate_regression.R"))
 #----------------------------------------------------------*
 # Multiple multivariate regression of genes vs weight loss in the field infections
 # Requires: hm, field
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, 
                         "analysis_multiple_multivariate_regression_field.R"))
 #----------------------------------------------------------*
 # Let's create a PCA showing the overlap of the field and laboratory
 # immune gene expression
 if (0) source(file.path(canalysis, 
                         "analysis_compare_PCA_lab_field.R"))
 #----------------------------------------------------------*
 # Panel 2: multivariate regression of genes vs weight loss in the field infections
 # Requires: hm, field
 #----------------------------------------------------------*
 if (0) source(file = "code/figure_creation/Panel_2.R")
 
 # ***********************************************************
 # Part 7: Using a random forest model to predict weight loss and extrapolate
 # to the wild
 # ***********************************************************
 #----------------------------------------------------------*
 # Creates random forest model: WL_predict_gene.RData
 #----------------------------------------------------------*
 if (1) source(file.path(canalysis, "analysis_random_forest_training.R"))
 #7.2: Aplication of random forest on field samples
 # requires hm and random forest model
 # Creates Field with updated new variable of predicted weight loss for each mouse
 # removed one mouse due to missing genotyping
 if (1) source(file.path(canalysis, "analysis_apply_random_forest.R"))
 # 7.3: We want to analyze the distribution of the predicted outcome variable
 # "WL_max" (predicted weight loss dependent on the immune gene expression values)
 # the distribution type is required for the downstream analysis
 # It seems that our predicted weight loss variale fits a normal distribution
 if (0) source(file.path(canalysis, "analysis_fit_distribution.R"))
 # 7.4: Here we test the hybrid impact on predicted weight loss and further test 
 # requires script: "analysis_apply_random_forest.R"
 # the effect of sex and infection status with Eimeria 
 # I have further tested to see if the presence of various parasites independent 
 # of eimeria have an impact on the predicted weight loss in combination with 
 # hybridicity 
 # # requires script: "analysis_apply_random_forest.R"
 if (0) source(file.path(canalysis, "analysis_tolerance.R"))
 # 7.7: Testing impact of infection status with Eimeria spp. on wild mice
 # Do we predict higher weight loss for mice that are infected with Eimeria spp? 
 # Is there any evidence for an association between infection status, 
 #hybridicity and their impact on weight loss? 
 if (0) source(file.path(canalysis, 
                         "analysis_predicted_WL_infection_with_Eimeria.R"))
 # 7.7: Testing impact of infection status with Eimeria spp. on wild mice
 # Do we predict higher weight loss for mice that are infected with Eimeria spp? 
 # Is there any evidence for an association between infection status, 
 #hybridicity and their impact on weight loss? 
 if (0) source(file.path(canalysis, 
                         "analysis_WL_parasites.R"))
 #----------------------------------------------------------*
 # Panel 3: Random Forest predicting weight loss
 # Requires: hm, field
 #----------------------------------------------------------*
 if (0) source(file = "code/figure_creation/Panel_3.R")
 #----------------------------------------------------------*
 #----------------------------------------------------------*
 # 7.8: Testing if CXCL9 is indeed important in the wild
 if (0) source(file.path(canalysis, 
                         "analysis_CXCL9.R"))
 
 
 
 

 
 
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
# Requires: previous script!!!
#----------------------------------------------------------*
if (0) source(file.path(canalysis, "heatmap_lab_genes.R"))
 #----------------------------------------------------------*
 # 6.3: Multiple multivariate regression of genes vs weight loss in the lab
 # Requires: hm, lab
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, "analysis_multiple_multivariate_regression.R"))
 #----------------------------------------------------------*
 # 6.4: Multiple multivariate regression of genes vs weight loss in the field infections
 # Requires: hm, field
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, "analysis_multiple_multivariate_regression_field.R"))
 
 
 # ***********************************************************
 # Part 7: Analysis                           ----
 # ***********************************************************
 #----------------------------------------------------------*
 #7.1: Random forest 
 # Training and testing a random forest that predicts weight loss
 # on experimental infections with Eimeria spp. 
 # Requires: hm
 # Creates random forest model: WL_predict_gene.RData
 #----------------------------------------------------------*
 if (0) source(file.path(canalysis, "analysis_random_forest_training.R"))
 #7.2: Aplication of random forest on field samples
 # requires hm and random forest model
 # Creates Field with updated new variable of predicted weight loss for each mouse
 # removed one mouse due to missing genotyping
 if (0) source(file.path(canalysis, "analysis_apply_random_forest.R"))
 # 7.3: We want to analyze the distribution of the predicted outcome variable
 # "WL_max" (predicted weight loss dependent on the immune gene expression values)
 # the distribution type is required for the downstream analysis
 # It seems that our predicted weight loss variale fits a normal distribution
 if (0) source(file.path(canalysis, "analysis_fit_distribution.R"))
 # 7.4: Here we test the hybrid impact on predicted weight loss and further test 
 # requires script: "analysis_apply_random_forest.R"
 # the effect of sex and infection status with Eimeria 
 # I have further tested to see if the presence of various parasites independent 
 # of eimeria have an impact on the predicted weight loss in combination with 
 # hybridicity 
 # # requires script: "analysis_apply_random_forest.R"
 if (0) source(file.path(canalysis, "analysis_hybrid_impact_sex_infection.R"))
 # 7.5: Exploring the relationships between variables and the predicted weight 
 # loss of field mice
 # Is there a relationship with the infection intensities with Eimeria spp.?
 # Is there a relationship between bmi and predicted weight loss?
 # We find a weak relationship of bmi and predicted weight loss
 if (0) source(file.path(canalysis, "analysis_bmi_predicted_WL.R"))
 # 7.6: Could we derive tolerance out of the predicted health impact and 
 # infection intenstities for each mouse? 
 # If so, can we test the hybridicity on this derives tolerance variable?
 # No evidence I can do that
 if (0) source(file.path(canalysis, "analysis_tolerance.R"))
 # 7.7: Testing impact of infection status with Eimeria spp. on wild mice
 # Do we predict higher weight loss for mice that are infected with Eimeria spp? 
 # Is there any evidence for an association between infection status, 
 #hybridicity and their impact on weight loss? 
 if (0) source(file.path(canalysis, 
                         "analysis_predicted_WL_infection_with_Eimeria.R"))
 # 7.7: Testing impact of infection status with Eimeria spp. on wild mice
 # Do we predict higher weight loss for mice that are infected with Eimeria spp? 
 # Is there any evidence for an association between infection status, 
 #hybridicity and their impact on weight loss? 
 if (0) source(file.path(canalysis, 
                         "analysis_WL_parasites.R"))
 ######
 if (0) source(file.path(canalysis, 
                         "analysis_multiple_multivariate_regression_field.R"))

 #analysis help next
 
# ***********************************************************
# Part end: Save information key information about project    ----
# ***********************************************************
sessionInfo()

 