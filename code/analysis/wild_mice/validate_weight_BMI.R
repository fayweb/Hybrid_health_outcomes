# =============================================================================
# Weight Loss and BMI Validation Analysis
# Purpose: Validate predicted weight loss using body weight and BMI metrics
# =============================================================================

# =============================================================================
# 1. BMI CALCULATION AND ANALYSIS
# =============================================================================

# Calculate BMI (Body Mass Index)
Field <- Field %>%
    dplyr::mutate(BMI = Body_Weight / (Body_Length / 100)^2)

# Check BMI distribution
shapiro.test(Field$BMI)  # Test for normality

# =============================================================================
# 2. BMI vs PREDICTED WEIGHT LOSS RELATIONSHIP
# =============================================================================

# Create scatter plot: BMI vs Predicted Weight Loss
bmi_plot <- ggplot(data = Field, aes(x = BMI, y = predicted_weight_loss)) +
    geom_point(color = "black", size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
    labs(
        title = "Body Mass Index vs Predicted Weight Loss",
        x = "Body Mass Index (BMI)",
        y = "Predicted Maximum Weight Loss (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
    )

# Display plot
print(bmi_plot)

# Fit linear model: BMI predicting weight loss
bmi_model <- lm(predicted_weight_loss ~ BMI, data = Field)
summary(bmi_model)
confint(bmi_model)

# Test correlation between BMI and predicted weight loss
cor_bmi <- cor.test(Field$BMI, Field$predicted_weight_loss, 
                    method = "spearman", use = "complete.obs")
print(cor_bmi)

# Generate model predictions for visualization
bmi_predictions <- ggpredict(bmi_model, terms = c("BMI"))

# Plot model predictions with confidence intervals
bmi_pred_plot <- ggplot(bmi_predictions, aes(x = x, y = predicted)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "blue") +
    labs(
        title = "Predicted Weight Loss based on BMI",
        x = "Body Mass Index (BMI)",
        y = "Predicted Weight Loss (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold")
    )

print(bmi_pred_plot)

# =============================================================================
# 3. BODY WEIGHT vs PREDICTED WEIGHT LOSS RELATIONSHIP
# =============================================================================

# Create scatter plot: Body Weight vs Predicted Weight Loss
weight_plot <- ggplot(data = Field, aes(x = Body_Weight, y = predicted_weight_loss)) +
    geom_point(color = "black", size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", size = 1) +
    labs(
        title = "Body Weight vs Predicted Weight Loss",
        x = "Body Weight (g)",
        y = "Predicted Maximum Weight Loss (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
    )

print(weight_plot)

# Fit linear model: Body Weight predicting weight loss
weight_model <- lm(predicted_weight_loss ~ Body_Weight, data = Field)
summary(weight_model)
confint(weight_model)

# Test correlation between body weight and predicted weight loss
cor_weight <- cor.test(Field$Body_Weight, Field$predicted_weight_loss, 
                       method = "spearman", use = "complete.obs")
print(cor_weight)

# =============================================================================
# 4. REVERSE RELATIONSHIP: PREDICTED WEIGHT LOSS vs BODY WEIGHT
# =============================================================================

# Model: How predicted weight loss relates to observed body weight
reverse_weight_model <- lm(Body_Weight ~ predicted_weight_loss, data = Field)
summary(reverse_weight_model)
confint(reverse_weight_model)

# Create scatter plot for reverse relationship
reverse_weight_plot <- ggplot(data = Field, aes(x = predicted_weight_loss, y = Body_Weight)) +
    geom_point(color = "black", size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgreen", size = 1) +
    labs(
        title = "Predicted Weight Loss vs Observed Body Weight",
        x = "Predicted Maximum Weight Loss (%)",
        y = "Body Weight (g)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
    )

print(reverse_weight_plot)

# =============================================================================
# 5. SUMMARY STATISTICS AND MODEL COMPARISON
# =============================================================================

# Create summary table of key statistics
validation_summary <- data.frame(
    Metric = c("BMI vs Predicted WL", "Body Weight vs Predicted WL", "Predicted WL vs Body Weight"),
    Correlation = c(cor_bmi$estimate, cor_weight$estimate, 
                    cor.test(Field$predicted_weight_loss, Field$Body_Weight, 
                             method = "spearman", use = "complete.obs")$estimate),
    P_value = c(cor_bmi$p.value, cor_weight$p.value,
                cor.test(Field$predicted_weight_loss, Field$Body_Weight, 
                         method = "spearman", use = "complete.obs")$p.value),
    R_squared = c(summary(bmi_model)$r.squared, 
                  summary(weight_model)$r.squared,
                  summary(reverse_weight_model)$r.squared)
)

print("Validation Summary:")
print(validation_summary)

# =============================================================================
# 6. SAVE PLOTS
# =============================================================================

# Save all plots
ggsave(filename = paste0(an_fi, "/BMI_vs_predicted_WL.jpeg"),
       plot = bmi_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/BMI_predictions.jpeg"),
       plot = bmi_pred_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/body_weight_vs_predicted_WL.jpeg"),
       plot = weight_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/predicted_WL_vs_body_weight.jpeg"),
       plot = reverse_weight_plot, width = 6, height = 4, dpi = 300)

# Save all plots
ggsave(filename = paste0(an_fi, "/BMI_vs_predicted_WL.pdf"),
       plot = bmi_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/BMI_predictions.pdf"),
       plot = bmi_pred_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/body_weight_vs_predicted_WL.pdf"),
       plot = weight_plot, width = 6, height = 4, dpi = 300)

ggsave(filename = paste0(an_fi, "/predicted_WL_vs_body_weight.pdf"),
       plot = reverse_weight_plot, width = 6, height = 4, dpi = 300)

# =============================================================================
# 7. MODEL DIAGNOSTICS (OPTIONAL)
# =============================================================================

# Check model assumptions for key models
par(mfrow = c(2, 2))

# BMI model diagnostics
plot(bmi_model, main = "BMI Model Diagnostics")

# Body weight model diagnostics  
plot(weight_model, main = "Body Weight Model Diagnostics")

# Reverse model diagnostics
plot(reverse_weight_model, main = "Reverse Weight Model Diagnostics")

par(mfrow = c(1, 1))  # Reset plotting parameters

cat("Analysis complete. Check plots and model summaries for validation results.\n")

ggsave(filename = paste0(an_fi, "/enhanced_weight_validation_v2.jpeg"),
       plot = enhanced_weight_plot_v2, width = 8, height = 6, dpi = 300)

