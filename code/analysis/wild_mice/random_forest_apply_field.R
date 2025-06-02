# ***********************************************************
# Script 7.1: Apply Random Forest Model to Wild Mice
# ***********************************************************
# Purpose: Apply laboratory-trained RF model to predict weight loss in wild mice
# Requires: Trained RF model, hm dataset with immune gene expression data
# Creates: hm dataset with predicted weight loss values, application summary stats
cat(strrep("=", 60), "\n")
cat("APPLYING RANDOM FOREST MODEL TO WILD MICE\n")
cat(strrep("=", 60), "\n")



# ---- Import Data ----
cat("Loading datasets...\n")

# Load the complete dataset
if (!exists("hm")) {
    hm <- read.csv("data/analysis/final/hm_ready_for_analysis.csv")
}

# Filter for field mice only
Field <- hm %>%
    filter(origin == "Field")

cat("✓ Loaded", nrow(Field), "wild-caught mice\n")

# ---- Load Random Forest Model ----
cat("Loading trained random forest model...\n")

# Load the trained model
weight_loss_predict <- readRDS("code/models/WL_predict_gene.rds")

cat("✓ Random forest model loaded successfully\n")
cat("  - Model type:", class(weight_loss_predict)[1], "\n")
cat("  - Number of trees:", weight_loss_predict$ntree, "\n")
cat("  - Variables used:", length(weight_loss_predict$forest$xlevels), "\n")

# ---- Prepare Gene Expression Data ----
cat("Preparing immune gene expression data for prediction...\n")

# Select gene columns for prediction
genes_for_prediction <- Field %>%
    ungroup() %>%
    dplyr::select(all_of(Genes_v))

cat("✓ Extracted", ncol(genes_for_prediction), "immune genes for", nrow(genes_for_prediction), "mice\n")

# Check for missing values
missing_genes <- colSums(is.na(genes_for_prediction))
if (any(missing_genes > 0)) {
    cat("⚠ Warning: Missing values detected in gene expression data:\n")
    print(missing_genes[missing_genes > 0])
} else {
    cat("✓ No missing values in gene expression data\n")
}

# ---- Apply Random Forest Model ----
cat("Generating weight loss predictions...\n")

# Set seed for reproducibility
set.seed(540)

# Generate predictions using the trained RF model
predicted_WL <- predict(weight_loss_predict, genes_for_prediction)

cat("✓ Predictions generated for", length(predicted_WL), "mice\n")

# ---- Add Predictions to Dataset ----
cat("Adding predictions to field dataset...\n")

# Add predictions to Field dataset
Field$predicted_weight_loss <- predicted_WL

# Also add to the main hm dataset for other scripts
hm$predicted_weight_loss <- NA
hm$predicted_weight_loss[hm$origin == "Field"] <- predicted_WL

cat("✓ Predictions added to datasets\n")

# ---- Summary Statistics ----
cat("PREDICTION SUMMARY STATISTICS:\n")
cat(strrep("-", 30), "\n")

prediction_summary <- summary(predicted_WL)
print(prediction_summary)

cat("\nPrediction distribution:\n")
cat("Mean ± SD:", round(mean(predicted_WL), 2), "±", round(sd(predicted_WL), 2), "\n")
cat("Range:", round(min(predicted_WL), 2), "to", round(max(predicted_WL), 2), "\n")
cat("IQR:", round(quantile(predicted_WL, 0.25), 2), "to", round(quantile(predicted_WL, 0.75), 2), "\n")

# Check for extreme predictions
extreme_low <- sum(predicted_WL < 0)
extreme_high <- sum(predicted_WL > 30)
cat("Extreme predictions: <0%:", extreme_low, " | >30%:", extreme_high, "\n")

# Continue with validation immediately
cat("\n", strrep("=", 60), "\n")
cat("STARTING FIELD VALIDATION ANALYSES\n")
cat(strrep("=", 60), "\n")

