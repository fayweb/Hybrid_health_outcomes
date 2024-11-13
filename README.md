# Hybrid_health_outcomes

Hybrid Health Outcomes
Overview
This repository, Hybrid_health_outcomes, contains the code and data for studying the health impacts of Eimeria infections on hybrid mice through predictive modeling and immune gene expression analysis. The goal is to translate laboratory findings into field applications, using immune markers to predict infection severity and health outcomes in wild populations.

Repository Structure
/code/ - Code files are organized by lab and field data processing, analysis, and model building:

lab: Scripts for data import, cleaning, and visualization of laboratory infection data.

field: Scripts for data import, processing, and analysis of field-collected infection data.

analysis: Core analyses, including PCA, linear regression, random forest training, and other statistical evaluations.

design: Contains experimental design details.

models: Model-building files, including random forest prediction scripts.

nmi: Data merging, normalization, and imputation scripts for combining lab and field data.

/data/ - Data files (organized by lab, field, and analysis datasets):

lab: Contains raw, intermediate, and final laboratory data.

field: Raw, intermediate, and processed field data.

analysis: Data products prepared for analysis.

/output/ - Output directory for figures, tables, and processed results:

figures: Various visualizations from the project, including gene expression and PCA plots.

tables: Statistical tables and regression outputs.

Master Script
The main script, Hybrid_health_outcomes.Rproj, controls the execution of all analysis steps, from data cleaning to predictive modeling. This script runs the following sections in sequence:

Environment Setup: Loads necessary packages and sets paths for project directories.

Data Cleaning: Imports raw data, formats it, and performs initial visualizations for lab and field data.

Data Merging, Normalization, and Imputation: Combines lab and field data, normalizes immune gene expression, and imputes missing values.

Experimental Design: Outlines key elements of the study, including sample distributions and parasite information.

Analysis:
PCA Analysis: Examines immune gene expression in lab data.

Linear Regressions: Analyzes relationships between gene expression components and health outcomes.

Random Forest: Trains and tests a model to predict weight loss based on immune markers.

Field Prediction: Applies the trained model to field data, assessing weight loss predictions in wild mice.

Key Scripts
functions.R - Custom functions for data processing, cleaning, and analysis.

analysis_PCA_genes_lab.R - Principal Component Analysis (PCA) on lab immune gene data.

analysis_random_forest_training.R - Random forest training and validation on lab infection data.

analysis_apply_random_forest.R - Application of random forest model to field data for health outcome predictions.

Getting Started

Clone the repository:

bash
Copy code
git clone https://github.com/fayweb/Hybrid_health_outcomes.git
cd Hybrid_health_outcomes
Install required R packages using pacman (included in setup script).

Run the Analysis: Use the main script (Hybrid_health_outcomes.Rproj) to execute all analyses or explore individual scripts in /code/ for specific tasks.

Results
The analysis highlights:

Predictive Model Accuracy: Random forest explains ~47% of the weight loss variance, with CXCL9 and TNF as significant immune markers.
Field Validation: Model predictions for wild mice align well with infection severity, supporting real-world application.
