# ===================================================================
# LINEAR REGRESSION ANALYSIS: PC AXES AS PREDICTORS OF WEIGHT LOSS
# ===================================================================
# Purpose: Test how immune gene expression patterns (PC1, PC2) predict
#          infection-induced weight loss in laboratory mice
# Author:Fay
# ==================================================================
# ===================================================================
# LINEAR REGRESSION ANALYSIS: PC AXES AS PREDICTORS OF WEIGHT LOSS
# ===================================================================
# Purpose: Test how immune gene expression patterns (PC1, PC2) predict
#          infection-induced weight loss in laboratory mice
# Author:Fay
# ==================================================================

# -------------------------------------------------------------------
# SECTION 1: CORE REGRESSION MODELS
# -------------------------------------------------------------------

# Model 1: PC axes only (simple immune signature model)
model_pc_only <- lm(WL_max ~ PC1 + PC2, data = lab)
cat("Model 1: PC axes only\n")
summary(model_pc_only)

# Model 2: Complete model with all covariates
model_complete <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE +
                         mouse_strain + immunization + weight_dpi0, data = lab)
cat("\nModel 2: Complete model\n")
summary(model_complete)

# Model 3: Interaction model (main biological interest)
model_interaction <- lm(WL_max ~ PC1 * current_infection + PC2 * current_infection, 
                        data = lab)
cat("\nModel 3: Interaction model\n")
summary(model_interaction)

# -------------------------------------------------------------------
# SECTION 2: MODEL COMPARISON TABLE
# -------------------------------------------------------------------

# Create model list
models_comparison <- list(
    "PC Only" = model_pc_only,
    "Complete" = model_complete, 
    "Interaction" = model_interaction
)

# Create clean comparison table
comparison_table <- modelsummary(
    models_comparison,
    output = "gt",
    stars = c('*' = .05, '**' = .01, '***' = .001),
    coef_map = c(
        "(Intercept)" = "Intercept",
        "PC1" = "PC1 (Inflammatory)",
        "PC2" = "PC2 (Regulatory)", 
        "current_infectionE. ferrisi" = "<i>E. ferrisi</i> infection",
        "current_infectionE. falciformis" = "<i>E. falciformis</i> infection",
        "delta_ct_cewe_MminusE" = "Infection intensity",
        "weight_dpi0" = "Initial body weight",
        "PC1:current_infectionE. ferrisi" = "PC1 × <i>E. ferrisi</i>",
        "PC1:current_infectionE. falciformis" = "PC1 × <i>E. falciformis</i>", 
        "current_infectionE. ferrisi:PC2" = "PC2 × <i>E. ferrisi</i>",
        "current_infectionE. falciformis:PC2" = "PC2 × <i>E. falciformis</i>"
    ),
    gof_map = c("nobs", "r.squared", "adj.r.squared", "statistic", "p.value"),
    notes = c("Reference: Uninfected controls", 
              "Mouse strain effects omitted for clarity")
) %>%
    tab_header(
        title = "Linear regression models: Immune signatures predict weight loss",
        subtitle = "Comparison of three modeling approaches"
    ) %>%
    tab_footnote(
        footnote = "PC1: inflammatory genes; PC2: regulatory genes"
    ) %>%
    cols_label(
        "PC Only" = "Immune Only",
        "Complete" = "Full Model", 
        "Interaction" = "Interaction"
    )

# Save comparison table
save_table_all_formats(comparison_table, "pca_regression_comparison")

# -------------------------------------------------------------------
# SECTION 3: SIMPLE INTERACTION MODEL TABLE
# -------------------------------------------------------------------

# Create a focused table for just the interaction model
interaction_only_table <- modelsummary(
    list("Interaction Model" = model_interaction),
    output = "gt",
    stars = c('*' = .05, '**' = .01, '***' = .001),
    coef_map = c(
        "(Intercept)" = "Intercept",
        "PC1" = "PC1 (Inflammatory)",
        "PC2" = "PC2 (Regulatory)", 
        "current_infectionE. ferrisi" = "E. ferrisi infection",
        "current_infectionE. falciformis" = "E. falciformis infection",
        "PC1:current_infectionE. ferrisi" = "PC1 × E. ferrisi",
        "PC1:current_infectionE. falciformis" = "PC1 × E. falciformis", 
        "current_infectionE. ferrisi:PC2" = "PC2 × E. ferrisi",
        "current_infectionE. falciformis:PC2" = "PC2 × E. falciformis"
    ),
    gof_map = c("nobs", "r.squared", "adj.r.squared", "statistic", "p.value"),
    notes = "Reference group: Uninfected controls"
) %>%
    tab_header(
        title = "Interaction Model: Immune Signatures × Parasite Species",
        subtitle = "Principal components predicting weight loss by infection group"
    )

# Save interaction table
save_table_all_formats(interaction_only_table, "interaction_model_detailed")

# -------------------------------------------------------------------
# SECTION 4: COEFFICIENT VISUALIZATION
# -------------------------------------------------------------------

# Create coefficient plot comparing all three models
coef_comparison_plot <- plot_summs(
    model_pc_only, model_complete, model_interaction,
    colors = c("#1f77b4", "#ff7f0e", "#2ca02c"),
    point.size = 3,
    model.names = c("PC Only", "Complete", "Interaction"),
    coefs = c("PC1", "PC2", "current_infectionE. ferrisi", "current_infectionE. falciformis",
              "PC1:current_infectionE. ferrisi", "PC1:current_infectionE. falciformis",
              "current_infectionE. ferrisi:PC2", "current_infectionE. falciformis:PC2")
) +
    labs(
        # title = "Coefficient comparison Across Models",
        x = "Coefficient Estimate",
        y = "Model Terms"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12)
    )

# Save coefficient plot
ggsave(
    filename = paste0(an_fi, "/coefficient_comparison_pca_models.jpeg"),
    plot = coef_comparison_plot,
    width = 10, height = 6, dpi = 300
)

# -------------------------------------------------------------------
# SECTION 5: INTERACTION EFFECT VISUALIZATIONS
# -------------------------------------------------------------------

# PC1 × Infection interaction plot
pc1_interaction_plot <- ggpredict(model_interaction, terms = c("PC1", "current_infection")) %>%
    plot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(
        title = NULL,
        x = "PC1 (Inflammatory Genes)",
        y = "Predicted Weight Loss (%)",
        color = "Infection Group"
    ) +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()
    )

# PC2 × Infection interaction plot  
pc2_interaction_plot <- ggpredict(model_interaction, terms = c("PC2", "current_infection")) %>%
    plot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(
        title = NULL,
        x = "PC2 (Regulatory Genes)",
        y = "Predicted Weight Loss (%)",
        color = "Infection Group"
    ) +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()
    )

# Save interaction plots
ggsave(
    filename = paste0(an_fi, "/pc1_infection_interaction.png"),
    plot = pc1_interaction_plot,
    width = 8, height = 6, dpi = 300
)

ggsave(
    filename = paste0(an_fi, "/pc2_infection_interaction.png"),
    plot = pc2_interaction_plot,
    width = 8, height = 6, dpi = 300
)

# -------------------------------------------------------------------
# SECTION 6: MODEL DIAGNOSTICS (OPTIONAL)
# -------------------------------------------------------------------

# Residual diagnostics for interaction model
qq_plot <- ggplot(data.frame(residuals = resid(model_interaction)), aes(sample = residuals)) +
    stat_qq(color = "steelblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(
        title = "Q-Q Plot: Residual Normality",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
    ) +
    theme_minimal()

residuals_plot <- ggplot(data.frame(
    fitted = fitted(model_interaction),
    residuals = resid(model_interaction)
), aes(x = fitted, y = residuals)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", color = "orange", se = FALSE) +
    labs(
        title = "Residuals vs Fitted Values",
        x = "Fitted Values",
        y = "Residuals"
    ) +
    theme_minimal()

# Save diagnostic plots
ggsave(
    filename = paste0(an_fi, "/model_diagnostics_qq.png"),
    plot = qq_plot,
    width = 6, height = 5, dpi = 300
)

ggsave(
    filename = paste0(an_fi, "/model_diagnostics_residuals.png"),
    plot = residuals_plot,
    width = 6, height = 5, dpi = 300
)

# -------------------------------------------------------------------
# SECTION 7: COMBINED PANEL FIGURE
# -------------------------------------------------------------------

# Create main results panel combining key visualizations
results_panel <- (pc1_interaction_plot | pc2_interaction_plot) /
    coef_comparison_plot +
    plot_layout(
        heights = c(2, 1),
        guides = 'collect'
    ) +
    plot_annotation(
        title = 'Principal Component Analysis: Immune Signatures Predict Weight Loss',
        tag_levels = 'A',
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )

# Save panel figure
ggsave(
    filename = paste0(panels_fi, "/pca_regression_results_panel.png"),
    plot = results_panel,
    width = 14, height = 10, dpi = 300
)

# -------------------------------------------------------------------
# SECTION 8: EXTRACT STATISTICS FOR TEXT
# -------------------------------------------------------------------

# Extract key statistics for reporting in text
extract_model_stats <- function(model) {
    s <- summary(model)
    list(
        r_squared = round(s$r.squared, 3),
        adj_r_squared = round(s$adj.r.squared, 3),
        f_stat = round(s$fstatistic[1], 2),
        p_value = round(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE), 6),
        n_obs = nobs(model)
    )
}

# Get stats for all models
pc_only_stats <- extract_model_stats(model_pc_only)
complete_stats <- extract_model_stats(model_complete)
interaction_stats <- extract_model_stats(model_interaction)

# Fix the statistics printing section:
cat("\n===============================================================\n")
cat("KEY STATISTICS FOR RESULTS TEXT\n")
cat("===============================================================\n")

cat(sprintf("PC Only Model: R² = %.3f, F = %.2f, p = %.6f, n = %d\n", 
            pc_only_stats$r_squared, pc_only_stats$f_stat, 
            pc_only_stats$p_value, pc_only_stats$n_obs))

cat(sprintf("Complete Model: R² = %.3f, F = %.2f, p = %.6f, n = %d\n", 
            complete_stats$r_squared, complete_stats$f_stat, 
            complete_stats$p_value, complete_stats$n_obs))

cat(sprintf("Interaction Model: R² = %.3f, F = %.2f, p = %.6f, n = %d\n", 
            interaction_stats$r_squared, interaction_stats$f_stat, 
            interaction_stats$p_value, interaction_stats$n_obs))

cat("===============================================================\n")
# -------------------------------------------------------------------
# SECTION 9: CLEAN UP
# -------------------------------------------------------------------

cat("\n✅ Analysis complete! Check the following folders:\n")
cat("📊 Tables:", file.path(tables), "\n")
cat("📈 Figures:", file.path(an_fi), "\n") 
cat("🎨 Panels:", file.path(panels_fi), "\n")

# Clean up intermediate objects (optional)
# rm(qq_plot, residuals_plot, pc_only_stats, complete_stats, interaction_stats)
