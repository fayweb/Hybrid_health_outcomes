# ***********************************************************
# Random Forest Model Diagnostics and Comparison
# Add this to your existing random forest script
# ***********************************************************
# ***********************************************************
# Part 1: Model Diagnostics for Test Set Model (43% model)
# ***********************************************************

# 1.1: Residuals vs Fitted Plot
residuals_data <- data.frame(
    fitted = result$predictions,
    residuals = result$WL_max - result$predictions,
    observed = result$WL_max
)

# Residuals vs Fitted
residuals_fitted_plot <- ggplot(residuals_data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.3) +
    labs(x = "Fitted Values", 
         y = "Residuals",
         title = "Residuals vs Fitted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# 1.2: Q-Q Plot of Residuals
qq_plot <- ggplot(residuals_data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(x = "Theoretical Quantiles", 
         y = "Sample Quantiles",
         title = "Q-Q Plot of Residuals") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# 1.3: Cross-validation Performance Curves
# Extract CV results from your existing predict_WL_cv object
cv_results <- data.frame(
    iteration = 1:length(predict_WL_cv$cv.stats$rmse),
    rmse = predict_WL_cv$cv.stats$rmse,
    var_exp = predict_WL_cv$cv.stats$var.exp,
    mae = predict_WL_cv$cv.stats$mae
)

# RMSE across CV iterations
cv_rmse_plot <- ggplot(cv_results, aes(x = iteration, y = rmse)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = mean(cv_results$rmse), 
               linetype = "dashed", color = "red") +
    labs(x = "Cross-Validation Iteration", 
         y = "RMSE",
         title = "Cross-Validation RMSE") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Variance explained across CV iterations
cv_var_plot <- ggplot(cv_results, aes(x = iteration, y = var_exp)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = mean(cv_results$var_exp), 
               linetype = "dashed", color = "red") +
    labs(x = "Cross-Validation Iteration", 
         y = "Variance Explained (%)",
         title = "Cross-Validation Variance Explained") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Combine diagnostic plots
diagnostic_panel <- grid.arrange(
    residuals_fitted_plot, qq_plot,
    cv_rmse_plot, cv_var_plot,
    ncol = 2, nrow = 2
)

# Save diagnostic plots
ggsave(filename = paste0(an_fi, "/rf_model_diagnostics.pdf"),
       plot = diagnostic_panel, width = 12, height = 10, dpi = 300)

# ***********************************************************
# Part 2: Model Comparison Table
# ***********************************************************

# Create comparison table between test model and final model
model_comparison <- data.frame(
    Model = c("Test Set Model", "Final Model"),
    Dataset = c("Training subset (70%)", "Complete dataset"),
    n_samples = c(nrow(train.data), nrow(gene_W)),
    n_trees = c(26, 308),
    Variance_Explained = c(43.55, 47.6),
    MSE = c(43.23, 32.47),
    RMSE = c(sqrt(43.23), sqrt(32.47)),
    Test_Correlation = c(0.79, "Applied to wild mice"),
    Purpose = c("Validation & performance assessment", "Wild mouse prediction")
)

# Print the table
print("Random Forest Model Comparison:")
print(model_comparison)

# Save as formatted table
write.csv(model_comparison, 
          file = paste0(tables, "/rf_model_comparison.csv"), 
          row.names = FALSE)

# Create a formatted display table
library(kableExtra)
comparison_table <- model_comparison %>%
    kbl(caption = "Random Forest Model Comparison",
        digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(1, background = "#f0f0f0") %>%
    row_spec(2, background = "#e6f3ff")

# Save formatted table
save_kable(comparison_table, 
           file = paste0(tables, "/rf_model_comparison_formatted.html"))

# ***********************************************************
# Part 3: Top Predictor Correlation Matrix
# ***********************************************************

# Get top 5 predictors based on importance
top_predictors <- var_imp %>%
    top_n(5, Importance) %>%
    pull(Genes)

# Create correlation matrix for top predictors
top_pred_data <- gene_W %>%
    dplyr::select(all_of(top_predictors))

# Calculate correlation matrix
cor_matrix <- cor(top_pred_data, use = "complete.obs")

# Create correlation plot
pdf(paste0(an_fi, "/top_predictors_correlation.pdf"), width = 8, height = 6)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.8,
         title = "Correlation Matrix: Top 5 Predictors",
         mar = c(0,0,2,0))
dev.off()

# Alternative ggplot version
cor_data <- as.data.frame(as.table(cor_matrix))
names(cor_data) <- c("Gene1", "Gene2", "Correlation")

cor_heatmap <- ggplot(cor_data, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "", title = "Top 5 Predictors: Correlation Matrix") +
    geom_text(aes(label = round(Correlation, 2)), size = 3)

ggsave(filename = paste0(an_fi, "/top_predictors_correlation_ggplot.pdf"),
       plot = cor_heatmap, width = 8, height = 6, dpi = 300)

# ***********************************************************
# Part 4: Summary Statistics
# ***********************************************************

# Model performance summary
cat("\n=== RANDOM FOREST MODEL SUMMARY ===\n")
cat("Test Set Performance (n =", nrow(test.data), "):\n")
cat("- Correlation (r):", round(cor(result$WL_max, result$predictions), 3), "\n")
cat("- R-squared:", round(cor(result$WL_max, result$predictions)^2, 3), "\n")
cat("- RMSE:", round(sqrt(mean((result$WL_max - result$predictions)^2)), 3), "\n")

cat("\nCross-Validation Results:\n")
cat("- Mean RMSE:", round(mean(cv_results$rmse), 3), "\n")
cat("- Mean Variance Explained:", round(mean(cv_results$var_exp), 3), "%\n")

cat("\nFinal Model (n =", nrow(gene_W), "):\n")
cat("- Variance Explained:", WL_predict_gene$rsq[length(WL_predict_gene$rsq)] * 100, "%\n")
cat("- MSE:", round(WL_predict_gene$mse[length(WL_predict_gene$mse)], 3), "\n")

cat("\nTop 5 Predictors:\n")
top_5 <- var_imp %>% 
    arrange(desc(Importance)) %>%
    head(5)
for(i in 1:nrow(top_5)) {
    cat(paste0(i, ". ", top_5$Genes[i], " (", round(top_5$Importance[i], 1), ")\n"))
}

cat("=====================================\n")

# ***********************************************************
# Part 5: Model Justification Text
# ***********************************************************

# Generate text for methods/results
cat("\n=== TEXT FOR MANUSCRIPT ===\n")
cat("Model Development Strategy:\n")
cat("We employed a two-stage random forest approach. First, we split the laboratory dataset")
cat("into training (70%) and testing (30%) sets to assess predictive performance and avoid")
cat("overfitting. After confirming robust performance (r = 0.79), we retrained the model")
cat("on the complete laboratory dataset to maximize predictive power for application to")
cat("wild-caught mice. This approach balances validation rigor with practical utility.\n")
cat("============================\n")