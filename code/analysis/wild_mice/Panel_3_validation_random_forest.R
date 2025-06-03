# ============================================================================
# FIGURE 3 PANEL CREATION - RANDOM FOREST FIELD VALIDATION
# ============================================================================
# Purpose: Create Figure 3 panel showing random forest field validation results
# Requires: Completed validation models from validate_infection_intensity.R
# Author: [Your name]
# Date: [Current date]
# ============================================================================


# ============================================================================
# LOAD DATA AND MODELS
# ============================================================================
cat("Loading data and models...\n")

# Load the validation results (adjust path as needed)
load("random_forest_validation_results.RData")



# ============================================================================
# PANEL C: Field Validation - Infection Intensity Effects  
# ============================================================================
cat("Creating Panel C: Infection intensity effects...\n")

# Use the actual validation results
coef_data <- data.frame(
    Predictor = c("Infected with Eimeria spp.", "Infection intensity with\nEimeria spp. (caecal qPCR)"),
    Estimate = c(validation_results$model_results_data$Estimate[4], 
                 validation_results$model_results_data$Estimate[5]),
    SE = c(validation_results$model_results_data$SE[4], 
           validation_results$model_results_data$SE[5]),
    p_value = c(validation_results$model_results_data$p_value[4], 
                validation_results$model_results_data$p_value[5])
) %>%
    mutate(
        conf.low = Estimate - 1.96 * SE,
        conf.high = Estimate + 1.96 * SE,
        significant = p_value < 0.05
    )

panel_C <- ggplot(coef_data, aes(x = Estimate, y = Predictor, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#4ACAFF")) +
    labs(x = "Estimate", y = "Predictor") +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 11)
    )

# ============================================================================
# PANEL D: Species-Specific Effects
# ============================================================================
cat("Creating Panel D: Species-specific effects...\n")

# Use the species model from validation results
species_preds <- ggpredict(validation_results$species_model, terms = "species_Eimeria")

panel_D <- ggplot(species_preds, aes(x = x, y = predicted, color = x)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
    scale_color_manual(values = color_mapping_species) +
    labs(
        x = "",
        y = "Predicted Weight Loss"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 11)
    )

# ============================================================================
# PANEL E: Infection Status Distribution
# ============================================================================
cat("Creating Panel E: Infection status raincloud plot...\n")

# Prepare data for raincloud plot
Field_plot <- Field %>%
    filter(!is.na(infection_status) & !is.na(predicted_weight_loss)) %>%
    mutate(
        infection_label = ifelse(infection_status, 
                                 "Infected with\nEimeria spp.", 
                                 "Eimeria spp.\nuninfected")
    )

panel_E <- Field_plot %>%
    ggplot(aes(x = predicted_weight_loss, y = infection_label, fill = infection_label)) +
    ggdist::stat_halfeye(
        adjust = 0.5,
        width = 0.6,
        alpha = 0.7,
        .width = 0,
        justification = -0.2,
        point_colour = NA
    ) +
    geom_boxplot(
        width = 0.15,
        outlier.shape = NA
    ) +
    ggdist::stat_dots(
        side = "left",
        justification = 1.1,
        binwidth = 0.25,
        alpha = 0.5
    ) +
    scale_fill_manual(values = color_mapping_infection) +
    labs(
        x = "Predicted weight loss",
        y = ""
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 11)
    )

# ============================================================================
# COMBINE PANELS INTO FIGURE 3
# ============================================================================
cat("Combining panels into Figure 3...\n")

# Create the panel layout
figure_3 <- (panel_A | panel_B) / 
    (panel_C | panel_D) / 
    panel_E +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(
        title = 'Figure 3',
        theme = theme(plot.title = element_text(size = 14, hjust = 0))
    )

# Display the figure
print(figure_3)

# ============================================================================
# SAVE FIGURE
# ============================================================================
cat("Saving Figure 3...\n")

ggsave(
    plot = figure_3,
    filename = paste0(panels_fi, "/Figure_3_RF_field_validation.pdf"),
    width = 12,
    height = 14,
    dpi = 300
)

ggsave(
    plot = figure_3,
    filename = paste0(panels_fi, "/Figure_3_RF_field_validation.png"),
    width = 12,
    height = 14,
    dpi = 300
)

# Also save individual panels if needed
ggsave(panel_A, filename = paste0(an_fi, "/panel_A_lab_validation.pdf"), width = 6, height = 5, dpi = 300)
ggsave(panel_B, filename = paste0(an_fi, "/panel_B_variable_importance.pdf"), width = 6, height = 5, dpi = 300)
ggsave(panel_C, filename = paste0(an_fi, "/panel_C_intensity_effects.pdf"), width = 6, height = 5, dpi = 300)
ggsave(panel_D, filename = paste0(an_fi, "/panel_D_species_effects.pdf"), width = 6, height = 5, dpi = 300)
ggsave(panel_E, filename = paste0(an_fi, "/panel_E_infection_distribution.pdf"), width = 8, height = 4, dpi = 300)

cat("✅ Figure 3 and individual panels saved successfully!\n")
cat("📁 Main figure: ", paste0(panels_fi, "/Figure_3_RF_field_validation.pdf"), "\n")
cat("📁 Individual panels saved in: ", an_fi, "\n")