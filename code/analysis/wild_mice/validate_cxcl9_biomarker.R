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

# 4. Test on 30% of lab data
lab_test_predictions <- predict(cxcl9_train_model, newdata = lab_test)
lab_test_correlation <- cor.test(lab_test_predictions, lab_test$WL_max)

# 5. Train final model on complete laboratory dataset
cxcl9_final_model <- lm(WL_max ~ CXCL9, data = lab)
final_summary <- summary(cxcl9_final_model)

# 6. Apply to field data
field_cxcl9_predictions <- predict(cxcl9_final_model, newdata = Field)

# 7. Compare CXCL9 predictions vs Random Forest predictions in field
cxcl9_vs_rf_correlation <- cor.test(field_cxcl9_predictions, Field$predicted_weight_loss)

# Print results
cat("Lab full model: R² =", round(final_summary$r.squared, 3), 
    ", p =", format(final_summary$coefficients[2,4], scientific = TRUE), "\n")
cat("Lab train-test validation: r =", round(lab_test_correlation$estimate, 3),
    ", p =", format(lab_test_correlation$p.value, scientific = TRUE), "\n")
cat("CXCL9 vs RF in field: r =", round(cxcl9_vs_rf_correlation$estimate, 3),
    ", p =", format(cxcl9_vs_rf_correlation$p.value, scientific = TRUE), "\n")


# Create CXCL9 validation summary table
cxcl9_validation_table <- data.frame(
    Analysis = c("Laboratory full model", "Laboratory train-test validation", 
                 "Field cross-population correlation"),
    n = c(136, 41, 336),
    Statistic = c("R²", "Pearson r", "Pearson r"),
    Value = c(0.128, 0.492, 0.172),
    p_value = c("< 0.001", "0.001", "0.002"),
    CI_95 = c("—", "[0.22, 0.70]", "[0.07, 0.27]")
)

# Convert to gt and save
cxcl9_gt_table <- cxcl9_validation_table %>%
    gt() %>%
    tab_header(title = "")

save_table_all_formats(cxcl9_gt_table, "Supplementary_Table_S12_CXCL9_validation")
