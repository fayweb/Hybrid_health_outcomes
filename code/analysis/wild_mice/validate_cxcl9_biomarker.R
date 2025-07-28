# Enhanced CXCL9 Model Performance Evaluation
# =====================================================

# 1. Test if CXCL9 alone is a good predictor in lab data
cxcl9_lab_full <- lm(WL_max ~ CXCL9, data = lab)
summary(cxcl9_lab_full)

# 2. Train-test split on laboratory data (70/30)
set.seed(333)  # Same seed as RF for consistency
train_indices <- sample(nrow(lab), size = 0.7 * nrow(lab))
lab_train <- lab[train_indices, ]
lab_test <- lab[-train_indices, ]

# 3. Train CXCL9 model on 70% of lab data
cxcl9_train_model <- lm(WL_max ~ CXCL9, data = lab_train)
train_summary <- summary(cxcl9_train_model)

# 4. Test on 30% of lab data with comprehensive performance metrics
lab_test_predictions <- predict(cxcl9_train_model, newdata = lab_test)

# *** NEW: Performance metrics for test set ***
# Correlation analysis
lab_test_correlation <- cor.test(lab_test_predictions, lab_test$WL_max)

# Root Mean Squared Error (RMSE)
rmse_test <- sqrt(mean((lab_test_predictions - lab_test$WL_max)^2))

# Mean Absolute Error (MAE)
mae_test <- mean(abs(lab_test_predictions - lab_test$WL_max))

# R-squared for test set (correlation squared)
r_squared_test <- cor(lab_test_predictions, lab_test$WL_max)^2

# *** NEW: Model diagnostic checks ***
# Check residuals normality (Shapiro-Wilk test)
residuals_train <- residuals(cxcl9_train_model)
shapiro_test <- shapiro.test(residuals_train)

# Check for heteroscedasticity (Breusch-Pagan test)
library(lmtest)
bp_test <- bptest(cxcl9_train_model)

# *** NEW: Compare performance with Random Forest ***
# Calculate RF performance metrics on same test set for comparison
# (Assuming you have RF predictions stored somewhere)
# rf_test_predictions <- predict(rf_model, lab_test)
# rf_rmse_test <- sqrt(mean((rf_test_predictions - lab_test$WL_max)^2))
# rf_r_test <- cor(rf_test_predictions, lab_test$WL_max)

# 5. Train final model on complete laboratory dataset
cxcl9_final_model <- lm(WL_max ~ CXCL9, data = lab)
final_summary <- summary(cxcl9_final_model)

# *** NEW: Final model performance metrics ***
final_rmse <- sqrt(mean(residuals(cxcl9_final_model)^2))
final_mae <- mean(abs(residuals(cxcl9_final_model)))

# 6. Apply to field data
field_cxcl9_predictions <- predict(cxcl9_final_model, newdata = Field)

# 7. Compare CXCL9 predictions vs Random Forest predictions in field
cxcl9_vs_rf_correlation <- cor.test(field_cxcl9_predictions, Field$predicted_weight_loss)

# *** NEW: Comprehensive results reporting ***
cat("=== CXCL9 MODEL PERFORMANCE SUMMARY ===\n\n")

cat("1. Laboratory Full Model (n=", nrow(lab), "):\n")
cat("   R² =", round(final_summary$r.squared, 3), "\n")
cat("   RMSE =", round(final_rmse, 3), "\n")
cat("   p-value =", format(final_summary$coefficients[2,4], scientific = TRUE), "\n\n")

cat("2. Train-Test Validation (test n=", nrow(lab_test), "):\n")
cat("   Correlation: r =", round(lab_test_correlation$estimate, 3), 
    " (95% CI:", round(lab_test_correlation$conf.int[1], 3), "-", 
    round(lab_test_correlation$conf.int[2], 3), ")\n")
cat("   p-value =", format(lab_test_correlation$p.value, scientific = TRUE), "\n")
cat("   R² =", round(r_squared_test, 3), "\n")
cat("   RMSE =", round(rmse_test, 3), "\n")
cat("   MAE =", round(mae_test, 3), "\n\n")

cat("3. Model Assumptions:\n")
cat("   Residuals normality (Shapiro-Wilk): p =", format(shapiro_test$p.value, scientific = TRUE), "\n")
cat("   Homoscedasticity (Breusch-Pagan): p =", format(bp_test$p.value, scientific = TRUE), "\n\n")

cat("4. Field Application (n=", nrow(Field), "):\n")
cat("   CXCL9 vs RF correlation: r =", round(cxcl9_vs_rf_correlation$estimate, 3), 
    " (95% CI:", round(cxcl9_vs_rf_correlation$conf.int[1], 3), "-", 
    round(cxcl9_vs_rf_correlation$conf.int[2], 3), ")\n")
cat("   p-value =", format(cxcl9_vs_rf_correlation$p.value, scientific = TRUE), "\n")

# *** NEW: Enhanced validation table with performance metrics ***
cxcl9_validation_table <- data.frame(
    Analysis = c("Laboratory full model", 
                 "Laboratory train-test validation", 
                 "Laboratory train-test validation",
                 "Laboratory train-test validation",
                 "Field cross-population correlation"),
    Metric = c("R²", "Pearson r", "RMSE", "MAE", "Pearson r"),
    Value = c(round(final_summary$r.squared, 3),
              round(lab_test_correlation$estimate, 3),
              round(rmse_test, 3),
              round(mae_test, 3),
              round(cxcl9_vs_rf_correlation$estimate, 3)),
    p_value = c(format(final_summary$coefficients[2,4], scientific = TRUE),
                format(lab_test_correlation$p.value, scientific = TRUE),
                "—",
                "—", 
                format(cxcl9_vs_rf_correlation$p.value, scientific = TRUE)),
    CI_95 = c("—", 
              paste0("[", round(lab_test_correlation$conf.int[1], 3), ", ", 
                     round(lab_test_correlation$conf.int[2], 3), "]"),
              "—",
              "—",
              paste0("[", round(cxcl9_vs_rf_correlation$conf.int[1], 3), ", ", 
                     round(cxcl9_vs_rf_correlation$conf.int[2], 3), "]")),
    n = c(136, 41, 41, 41, 336)
)

# *** NEW: Diagnostic plots for methods validation ***
par(mfrow = c(2, 2))

# 1. Observed vs Predicted (test set)
plot(lab_test$WL_max, lab_test_predictions,
     xlab = "Observed Weight Loss (%)", 
     ylab = "CXCL9-Predicted Weight Loss (%)",
     main = "CXCL9 Model: Observed vs Predicted (Test Set)")
abline(0, 1, col = "red", lty = 2)
text(min(lab_test$WL_max), max(lab_test_predictions), 
     paste("r =", round(lab_test_correlation$estimate, 3)), adj = 0)

# 2. Residuals vs Fitted
plot(fitted(cxcl9_train_model), residuals(cxcl9_train_model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# 3. Q-Q plot for residuals
qqnorm(residuals(cxcl9_train_model), main = "Normal Q-Q Plot")
qqline(residuals(cxcl9_train_model), col = "red")

# 4. CXCL9 vs RF predictions in field
plot(Field$predicted_weight_loss, field_cxcl9_predictions,
     xlab = "Random Forest Predictions (%)",
     ylab = "CXCL9 Predictions (%)",
     main = "Field Data: CXCL9 vs Random Forest")
abline(0, 1, col = "blue", lty = 2)
text(min(Field$predicted_weight_loss), max(field_cxcl9_predictions),
     paste("r =", round(cxcl9_vs_rf_correlation$estimate, 3)), adj = 0)

par(mfrow = c(1, 1))

# Convert to gt and save enhanced table
library(gt)
# Improved CXCL9 validation table formatting
cxcl9_validation_table <- data.frame(
    Analysis = c("Laboratory full model", 
                 "Laboratory train-test validation", 
                 "Laboratory train-test validation",
                 "Laboratory train-test validation",
                 "Field cross-population correlation"),
    Metric = c("R²", "Pearson r", "RMSE", "MAE", "Pearson r"),
    Value = c("0.128", "0.492", "6.36", "5.29", "0.172"),
    p_value = c("< 0.001", "0.001", "—", "—", "0.002"),
    CI_95 = c("—", "[0.22, 0.70]", "—", "—", "[0.07, 0.27]"),
    n = c(136, 41, 41, 41, 336)
)

# Better gt table formatting
library(gt)
cxcl9_gt_table <- cxcl9_validation_table %>%
    gt() %>%
    tab_header(
        title = "CXCL9 Single-Gene Model Performance Validation"
    ) %>%
    cols_label(
        Analysis = "Analysis",
        Metric = "Metric", 
        Value = "Value",
        p_value = "p-value",
        CI_95 = "95% CI",
        n = "n"
    ) %>%
    tab_footnote(
        footnote = "RMSE = Root Mean Squared Error; MAE = Mean Absolute Error",
        locations = cells_column_labels(columns = Metric)
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    )

# Save the improved table
save_table_all_formats(cxcl9_gt_table, "Supplementary_Table_CXCL9_validation_final")
