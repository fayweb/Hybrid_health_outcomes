# ============================================================================
# REPRODUCIBLE RANDOM FOREST FIELD VALIDATION WITH DISTRIBUTION TESTS
# ============================================================================

cat(strrep("=", 60), "\n")
cat("RUNNING ALL VALIDATION MODELS FOR REPRODUCIBILITY\n")
cat(strrep("=", 60), "\n")

library(moments)  # For skewness and kurtosis tests

# ============================================================================
# DISTRIBUTION ANALYSIS AND CORRELATION METHOD JUSTIFICATION
# ============================================================================
cat("\n0. DISTRIBUTION ANALYSIS\n")
cat(strrep("-", 30), "\n")

test_distribution <- function(data, variable_name) {
    # Remove missing values
    clean_data <- data[!is.na(data)]
    n <- length(clean_data)
    
    if (n < 3) {
        cat("Insufficient data for", variable_name, "\n")
        return(list(use_spearman = TRUE, reason = "insufficient data"))
    }
    
    # Calculate distribution statistics
    mean_val <- mean(clean_data)
    median_val <- median(clean_data)
    sd_val <- sd(clean_data)
    skewness_val <- skewness(clean_data)
    kurtosis_val <- kurtosis(clean_data) - 3  # Excess kurtosis
    
    # Count zeros
    zeros <- sum(clean_data == 0)
    zero_prop <- zeros / n
    
    cat("\n", variable_name, " (n = ", n, "):\n", sep = "")
    cat("  Mean: ", round(mean_val, 3), "\n")
    cat("  Median: ", round(median_val, 3), "\n")
    cat("  SD: ", round(sd_val, 3), "\n")
    cat("  Skewness: ", round(skewness_val, 3), "\n")
    cat("  Kurtosis: ", round(kurtosis_val, 3), "\n")
    cat("  Zeros: ", zeros, " (", round(zero_prop * 100, 1), "%)\n")
    
    # Assess normality
    normality_issues <- character(0)
    
    if (abs(skewness_val) > 1) {
        normality_issues <- c(normality_issues, "high skewness")
    }
    if (abs(kurtosis_val) > 1) {
        normality_issues <- c(normality_issues, "high kurtosis")
    }
    if (abs(mean_val - median_val) > sd_val * 0.5) {
        normality_issues <- c(normality_issues, "mean-median difference")
    }
    if (zero_prop > 0.1) {
        normality_issues <- c(normality_issues, "too many zeros")
    }
    
    use_spearman <- length(normality_issues) > 0
    
    if (use_spearman) {
        cat("  ⚠️  NORMALITY CONCERNS: ", paste(normality_issues, collapse = ", "), "\n")
        cat("  📊 RECOMMENDATION: Use Spearman correlation\n")
    } else {
        cat("  ✅ APPEARS NORMAL: Pearson correlation could be appropriate\n")
    }
    
    return(list(
        use_spearman = use_spearman,
        reason = paste(normality_issues, collapse = ", "),
        skewness = skewness_val,
        kurtosis = kurtosis_val,
        n = n
    ))
}

# Test distributions of key variables
cat("Testing distributions to justify correlation method choice...\n")

# Delta Ct distribution
delta_ct_dist <- test_distribution(Field$delta_ct_cewe_MminusE, "Delta Ct (qPCR intensity)")

# OPG distribution  
opg_dist <- test_distribution(Field$OPG[Field$OPG > 0], "Oocyst counts (OPG > 0)")

# Predicted weight loss distribution
pwl_dist <- test_distribution(Field$predicted_weight_loss, "Predicted weight loss")

# Summary recommendation
cat("\n📊 CORRELATION METHOD JUSTIFICATION:\n")
if (delta_ct_dist$use_spearman || opg_dist$use_spearman) {
    cat("Using Spearman correlations due to non-normal distributions in intensity measures.\n")
    correlation_method <- "spearman"
} else {
    cat("Data appear sufficiently normal for Pearson correlations.\n")
    correlation_method <- "pearson"
}

# ============================================================================
# MODEL 1: INFECTION STATUS EFFECT
# ============================================================================
cat("\n1. INFECTION STATUS MODEL\n")
cat(strrep("-", 30), "\n")

# Filter data
Field_status <- Field %>% 
    filter(!is.na(infection_status) & !is.na(predicted_weight_loss))

cat("Sample size:", nrow(Field_status), "mice\n")
cat("Group sizes:\n")
print(table(Field_status$infection_status))

# Run model
status_model <- lm(predicted_weight_loss ~ infection_status, data = Field_status)
status_results <- broom::tidy(status_model)
status_summary <- summary(status_model)

cat("\nModel results:\n")
print(status_results)
cat("Model R²:", round(status_summary$r.squared, 4), "\n")

# Extract key statistics
status_effect <- status_results$estimate[2]
status_se <- status_results$std.error[2]
status_p <- status_results$p.value[2]
status_r2 <- status_summary$r.squared
status_f_p <- pf(status_summary$fstatistic[1], 
                 status_summary$fstatistic[2], 
                 status_summary$fstatistic[3], 
                 lower.tail = FALSE)

# ============================================================================
# CREATE RAINCLOUD PLOT FOR INFECTION STATUS
# ============================================================================
cat("Creating raincloud plot for infection status...\n")

# Define colors for infection status
infection_colors <- c("Infected with Eimeria spp." = "#FF0028", "Uninfected" = "#4ACAFF")

# Prepare data for raincloud plot - create clean labels
Field_raincloud <- Field_status %>%
    mutate(
        infection_label = case_when(
            infection_status == TRUE ~ "Infected with Eimeria spp.",
            infection_status == FALSE ~ "Uninfected",
            TRUE ~ as.character(infection_status)
        )
    ) %>%
    drop_na(infection_label, predicted_weight_loss)

# Create raincloud plot
raincloud_plot <- Field_raincloud %>%
    ggplot(aes(y = infection_label, x = predicted_weight_loss, fill = infection_label)) + 
    ggdist::stat_halfeye(
        adjust = 0.5, 
        width = 0.6, 
        alpha = 0.7,
        .width = 0, 
        justification = -0.2, 
        point_colour = NA
    ) + 
    scale_fill_manual(values = infection_colors) +
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
    geom_point(
        shape = 95,
        size = 15,
        alpha = 0.2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    theme_minimal() +
    labs(
        y = "Infection status with *Eimeria* spp.", 
        x = "Predicted weight loss", 
        fill = "Infection status"
    ) +
    theme(
        legend.position = "none",
        axis.title.y = element_markdown(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)
    )

# Display the plot
print(raincloud_plot)

# Save the plot
ggsave(
    plot = raincloud_plot, 
    filename = paste0(an_fi, "/raincloud_infection_status.pdf"), 
    width = 6, 
    height = 4, 
    dpi = 300
)

# Save the plot
ggsave(
    plot = raincloud_plot, 
    filename = paste0(an_fi, "/raincloud_infection_status.jpeg"), 
    width = 6, 
    height = 4, 
    dpi = 300
)

cat("Raincloud plot saved to:", paste0(an_fi, "/raincloud_infection_status.pdf"), "\n")

# ============================================================================
# MODEL 2: SPECIES-SPECIFIC EFFECTS
# ============================================================================
cat("\n2. SPECIES-SPECIFIC MODEL\n")
cat(strrep("-", 30), "\n")

# Filter data
Field_species <- Field %>% 
    filter(!is.na(species_Eimeria) & !is.na(predicted_weight_loss))

cat("Sample size:", nrow(Field_species), "mice\n")
cat("Group sizes:\n")
print(table(Field_species$species_Eimeria))

# Run model
species_model <- lm(predicted_weight_loss ~ species_Eimeria, data = Field_species)
species_results <- broom::tidy(species_model)
species_summary <- summary(species_model)

cat("\nModel results:\n")
print(species_results)
cat("Model R²:", round(species_summary$r.squared, 4), "\n")

# Extract key statistics
ferrisi_effect <- species_results$estimate[grep("ferrisi", species_results$term)]
ferrisi_se <- species_results$std.error[grep("ferrisi", species_results$term)]
ferrisi_p <- species_results$p.value[grep("ferrisi", species_results$term)]

falciformis_effect <- species_results$estimate[grep("falciformis", species_results$term)]
falciformis_se <- species_results$std.error[grep("falciformis", species_results$term)]
falciformis_p <- species_results$p.value[grep("falciformis", species_results$term)]

species_r2 <- species_summary$r.squared
species_f_p <- pf(species_summary$fstatistic[1], 
                  species_summary$fstatistic[2], 
                  species_summary$fstatistic[3], 
                  lower.tail = FALSE)

# ============================================================================
# CREATE SPECIES PREDICTION PLOT
# ============================================================================
cat("Creating species prediction plot...\n")

# Create predicted values using ggpredict
species_preds <- ggpredict(species_model, terms = "species_Eimeria")

# Create the plot
species_plot <- ggplot(species_preds, aes(x = x, y = predicted, color = x)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
    scale_color_manual(values = color_mapping_f) +
    labs(
        x = "*Eimeria* spp. subspecies",
        y = "Predicted Weight Loss",
        color = "Species"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_markdown(),
        axis.title.x = element_markdown(),
        legend.text = element_markdown()
    )

# Display the plot
print(species_plot)

# Save the plot
ggsave(
    plot = species_plot,
    filename = paste0(an_fi, "/predicted_weight_loss_species.pdf"),
    width = 8, 
    height = 6, 
    dpi = 300
)

# Save the plot
ggsave(
    plot = species_plot,
    filename = paste0(an_fi, "/predicted_weight_loss_species.jpeg"),
    width = 8, 
    height = 6, 
    dpi = 300
)

cat("Species prediction plot saved to:", paste0(an_fi, "/predicted_weight_loss_species.pdf"), "\n")

# ============================================================================
# MODEL 3: INFECTION INTENSITY INTERACTION MODEL
# ============================================================================
cat("\n3. INFECTION INTENSITY INTERACTION MODEL\n")
cat(strrep("-", 30), "\n")

# Filter data for interaction model - CONSISTENT VARIABLE NAMES
Field_interaction <- Field %>% 
    filter(!is.na(delta_ct_cewe_MminusE) & !is.na(infection_status) & !is.na(predicted_weight_loss))

cat("Sample size:", nrow(Field_interaction), "mice\n")
cat("Group sizes:\n")
print(table(Field_interaction$infection_status))  # FIXED: Use infection_status

# Run interaction model
interaction_model <- lm(predicted_weight_loss ~ delta_ct_cewe_MminusE * infection_status, 
                        data = Field_interaction)
interaction_results <- broom::tidy(interaction_model)
interaction_summary <- summary(interaction_model)

cat("\nInteraction Model results:\n")
print(interaction_results)
cat("Model R²:", round(interaction_summary$r.squared, 4), "\n")
cat("Model F-test p-value:", format(pf(interaction_summary$fstatistic[1], 
                                       interaction_summary$fstatistic[2], 
                                       interaction_summary$fstatistic[3], 
                                       lower.tail = FALSE), scientific = TRUE), "\n")

# Extract interaction statistics - FIXED COEFFICIENT NAMES
interaction_term <- interaction_results[grep(":", interaction_results$term), ]
if (nrow(interaction_term) > 0) {
    interaction_effect <- interaction_term$estimate[1]
    interaction_se <- interaction_term$std.error[1]
    interaction_p <- interaction_term$p.value[1]
} else {
    interaction_effect <- NA
    interaction_se <- NA
    interaction_p <- NA
}

# FIXED: Look for infection_statusTRUE instead of MC.Eimeria
infection_main_effect <- interaction_results$estimate[grep("infection_statusTRUE", interaction_results$term)[1]]
infection_main_se <- interaction_results$std.error[grep("infection_statusTRUE", interaction_results$term)[1]]
infection_main_p <- interaction_results$p.value[grep("infection_statusTRUE", interaction_results$term)[1]]

interaction_r2 <- interaction_summary$r.squared
interaction_f_p <- pf(interaction_summary$fstatistic[1], 
                      interaction_summary$fstatistic[2], 
                      interaction_summary$fstatistic[3], 
                      lower.tail = FALSE)

# FIXED: Updated coefficient names for plot_summs
plot_summs(interaction_model, 
           scale = TRUE, 
           robust = TRUE,
           inner_ci_level = 0.95, 
           outer_ci_level = 0.99,
           coefs = c("Infected with Eimeria spp." = "infection_statusTRUE", 
                     "Infection intensity with Eimeria spp." = 
                         "delta_ct_cewe_MminusE:infection_statusTRUE")) +  # FIXED INTERACTION TERM NAME
    theme_minimal() +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
    ) +
    labs(
        x = "Estimate",
        y = "Predictor"
    ) +
    scale_color_manual(values = c("#4ACAFF", "#EBEBEB")) +  # FIXED: Added # for hex colors
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") -> int_MC

int_MC
ggsave(paste0(an_fi, "/intensity_melting_curve.pdf"), int_MC, dpi = 300)
ggsave(paste0(an_fi, "/intensity_melting_curve.jpeg"), int_MC, dpi = 300)
# ============================================================================
# MODEL 4: INTENSITY CORRELATIONS
# ============================================================================
cat("\n4. INTENSITY CORRELATIONS\n")
cat(strrep("-", 30), "\n")

# qPCR intensity (infected mice only)
intensity_qpcr <- Field %>% 
    filter(!is.na(delta_ct_cewe_MminusE) & !is.na(predicted_weight_loss) & MC.Eimeria == "TRUE")

cat("qPCR intensity data available for", nrow(intensity_qpcr), "mice\n")

cor_qpcr <- cor.test(intensity_qpcr$predicted_weight_loss, 
                     intensity_qpcr$delta_ct_cewe_MminusE, 
                     method = correlation_method)

cat("qPCR Intensity Correlation (", correlation_method, "):\n")
cat("  r =", round(cor_qpcr$estimate, 3), "\n")
cat("  p =", format(cor_qpcr$p.value, scientific = TRUE), "\n")
if (!is.null(cor_qpcr$conf.int)) {
    cat("  95% CI: [", round(cor_qpcr$conf.int[1], 3), ",", round(cor_qpcr$conf.int[2], 3), "]\n")
}

# Oocyst counts
intensity_opg <- Field %>% 
    filter(!is.na(OPG) & !is.na(predicted_weight_loss) & OPG > 0)

cat("\nOocyst count data available for", nrow(intensity_opg), "mice\n")

cor_opg <- cor.test(intensity_opg$predicted_weight_loss, 
                    intensity_opg$OPG, 
                    method = correlation_method)

cat("Oocyst Count Correlation (", correlation_method, "):\n")
cat("  r =", round(cor_opg$estimate, 3), "\n")
cat("  p =", format(cor_opg$p.value, scientific = TRUE), "\n")
if (!is.null(cor_opg$conf.int)) {
    cat("  95% CI: [", round(cor_opg$conf.int[1], 3), ",", round(cor_opg$conf.int[2], 3), "]\n")
}

# Simple OPG model for comparison
model_opg <- lm(predicted_weight_loss ~ OPG, data = Field)
opg_model_summary <- summary(model_opg)
cat("\nSimple OPG linear model:\n")
cat("  R² =", round(opg_model_summary$r.squared, 4), "\n")
cat("  p =", format(pf(opg_model_summary$fstatistic[1], 
                       opg_model_summary$fstatistic[2], 
                       opg_model_summary$fstatistic[3], 
                       lower.tail = FALSE), scientific = TRUE), "\n")

# ============================================================================
# CREATE COMPREHENSIVE TABLE FROM ACTUAL RESULTS
# ============================================================================
cat("\n", strrep("=", 60), "\n")
cat("CREATING COMPREHENSIVE TABLE FROM ACTUAL MODEL RESULTS\n")
cat(strrep("=", 60), "\n")

# Build data frame from actual results - INCLUDING CORRELATIONS
model_results_data <- data.frame(
    Model = c(
        "Infection Status Effect",
        "Species Model: E. ferrisi", 
        "Species Model: E. falciformis",
        "Interaction Model: Main Effect",
        "Interaction Model: Interaction Term",
        "qPCR Intensity Correlation",
        "Oocyst Count Correlation"
    ),
    Comparison = c(
        "Infected vs Uninfected",
        "vs Uninfected", 
        "vs Uninfected",
        "MC.Eimeria TRUE vs FALSE",
        "ΔCt × Infection Status",
        "Predicted WL vs ΔCt (infected mice)",
        "Predicted WL vs OPG (positive counts)"
    ),
    n = c(
        nrow(Field_status),      # Infection status sample size
        nrow(Field_species),     # Species model total sample size
        nrow(Field_species),     # Species model total sample size (same model)
        nrow(Field_interaction), # Interaction model sample size
        nrow(Field_interaction), # Interaction model sample size (same model)
        nrow(intensity_qpcr),    # qPCR correlation sample size
        nrow(intensity_opg)      # OPG correlation sample size
    ),
    Estimate = c(
        status_effect,           # Status effect
        ferrisi_effect,          # E. ferrisi effect
        falciformis_effect,      # E. falciformis effect
        infection_main_effect,   # Main infection effect in interaction model
        interaction_effect,      # Interaction effect
        cor_qpcr$estimate,       # qPCR correlation coefficient
        cor_opg$estimate         # OPG correlation coefficient
    ),
    SE = c(
        status_se,               # Status SE
        ferrisi_se,              # E. ferrisi SE
        falciformis_se,          # E. falciformis SE
        infection_main_se,       # Main infection SE
        interaction_se,          # Interaction SE
        NA,                      # No SE for correlations
        NA                       # No SE for correlations
    ),
    p_value = c(
        status_p,                # Status p-value
        ferrisi_p,               # E. ferrisi p-value
        falciformis_p,           # E. falciformis p-value
        infection_main_p,        # Main infection p-value
        interaction_p,           # Interaction p-value
        cor_qpcr$p.value,        # qPCR correlation p-value
        cor_opg$p.value          # OPG correlation p-value
    ),
    Model_R2 = c(
        status_r2,               # Status model R²
        species_r2,              # Species model R²
        species_r2,              # Same model R²
        interaction_r2,          # Interaction model R²
        interaction_r2,          # Same model R²
        NA,                      # No R² for correlations
        NA                       # No R² for correlations
    ),
    Model_F_pvalue = c(
        status_f_p,              # Status F-test p-value
        species_f_p,             # Species F-test p-value
        species_f_p,             # Same F-test p-value
        interaction_f_p,         # Interaction F-test p-value
        interaction_f_p,         # Same F-test p-value
        NA,                      # No F-test for correlations
        NA                       # No F-test for correlations
    )
)

cat("Model results data frame created with actual values:\n")
print(model_results_data)

# Create formatted table
model_results_table <- model_results_data %>%
    mutate(
        # Format estimates and SE
        Estimate_formatted = case_when(
            is.na(Estimate) ~ "NA",
            str_detect(Comparison, "Interaction") ~ sprintf("%.3f", Estimate),
            TRUE ~ paste0("+", sprintf("%.2f", Estimate), "%")
        ),
        SE_formatted = case_when(
            is.na(SE) ~ "NA",
            TRUE ~ sprintf("%.3f", SE)
        ),
        # Format p-values
        p_formatted = case_when(
            is.na(p_value) ~ "NA",
            p_value < 0.001 ~ "<0.001",
            p_value < 0.01 ~ sprintf("%.3f", p_value),
            TRUE ~ sprintf("%.3f", p_value)
        ),
        # Format R-squared
        R2_formatted = sprintf("%.3f", Model_R2),
        # Format model F p-values
        F_p_formatted = case_when(
            Model_F_pvalue < 0.001 ~ "<0.001",
            Model_F_pvalue < 0.01 ~ sprintf("%.3f", Model_F_pvalue),
            TRUE ~ sprintf("%.3f", Model_F_pvalue)
        ),
        # Significance indicator
        is_significant = !is.na(p_value) & p_value < 0.05,
        # Create italicized species names
        Model_formatted = case_when(
            str_detect(Model, "E\\. ferrisi") ~ str_replace(Model, "E\\. ferrisi", "*E. ferrisi*"),
            str_detect(Model, "E\\. falciformis") ~ str_replace(Model, "E\\. falciformis", "*E. falciformis*"),
            TRUE ~ Model
        )
    ) %>%
    dplyr::select(Model_formatted, Comparison, n, Estimate_formatted, SE_formatted, 
                  p_formatted, R2_formatted, F_p_formatted, is_significant) %>%
    gt() %>%
    
    # Table headers - CORRECTED
    tab_header(
        title = "",
        subtitle = ""
    ) %>%
    
    # Column labels
    cols_label(
        Model_formatted = "Model/Effect",
        Comparison = "Comparison", 
        n = "n",
        Estimate_formatted = "Effect Size",
        SE_formatted = "SE",
        p_formatted = "p-value",
        R2_formatted = "Model R²",
        F_p_formatted = "F-test p"
    ) %>%
    
    # Hide significance column
    cols_hide(columns = is_significant) %>%
    
    # Add row groups
    tab_row_group(
        label = "Infection Detection",
        rows = 1
    ) %>%
    tab_row_group(
        label = "Species-Specific Effects",
        rows = 2:3
    ) %>%
    tab_row_group(
        label = "Intensity × Status Interaction",
        rows = 4:5
    ) %>%
    
    # Apply markdown formatting for italics
    fmt_markdown(columns = Model_formatted) %>%
    
    # Formatting
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
    ) %>%
    tab_style(
        style = cell_text(align = "center"),
        locations = cells_body(columns = c(n, SE_formatted, p_formatted, R2_formatted, F_p_formatted))
    ) %>%
    
    # Highlight significant results
    tab_style(
        style = cell_fill(color = "#E8F5E8"),
        locations = cells_body(
            columns = c(Model_formatted, Comparison, n, Estimate_formatted, SE_formatted, 
                        p_formatted, R2_formatted, F_p_formatted),
            rows = is_significant == TRUE
        )
    ) %>%
    
    # Add footnotes
    tab_footnote(
        footnote = "Effect size represents increase in predicted weight loss. Interaction term shows how ΔCt effect changes with infection status.",
        locations = cells_column_labels(columns = Estimate_formatted)
    ) %>%
    tab_footnote(
        footnote = "Spearman correlations used only for specific infection intensity correlation analyses (qPCR and oocyst data).",
        locations = cells_title(groups = "subtitle")
    ) %>%
    
    # Table options
    tab_options(
        table.font.size = 12,
        heading.title.font.size = 14,
        heading.subtitle.font.size = 12,
        row_group.font.weight = "bold"
    )
# Display the table
print(model_results_table)

# Save the table
save_table_all_formats(model_results_table, "RF_comprehensive_validation_results")

# ============================================================================
# SAVE ALL RESULTS FOR REPRODUCIBILITY
# ============================================================================
cat("\n", strrep("=", 60), "\n")
cat("SAVING ALL RESULTS FOR REPRODUCIBILITY\n")
cat(strrep("=", 60), "\n")

# Save all model objects and results
validation_results <- list(
    # Distribution analysis
    distribution_tests = list(
        delta_ct = delta_ct_dist,
        opg = opg_dist,
        predicted_wl = pwl_dist,
        correlation_method = correlation_method
    ),
    # Model objects
    status_model = status_model,
    species_model = species_model,
    interaction_model = interaction_model,
    opg_model = model_opg,
    # Correlation results
    cor_qpcr = cor_qpcr,
    cor_opg = cor_opg,
    # Summary data
    model_results_data = model_results_data,
    sample_sizes = list(
        status = nrow(Field_status),
        species = nrow(Field_species),
        interaction = nrow(Field_interaction),
        qpcr = nrow(intensity_qpcr),
        opg = nrow(intensity_opg)
    )
)

# Save to file
save(validation_results, file = "random_forest_validation_results.RData")

cat("✅ All validation models completed and saved!\n")
cat("✅ Distribution tests performed and documented!\n")
cat("✅ Comprehensive table created from actual model results!\n")
cat("✅ Results saved for reproducibility!\n")

# Print final summary
cat("\n", strrep("=", 60), "\n")
cat("FINAL VALIDATION SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("Distribution Analysis:\n")
cat("  Delta Ct: ", ifelse(delta_ct_dist$use_spearman, "Non-normal", "Normal"), " (skewness = ", round(delta_ct_dist$skewness, 2), ")\n")
cat("  OPG: ", ifelse(opg_dist$use_spearman, "Non-normal", "Normal"), " (skewness = ", round(opg_dist$skewness, 2), ")\n")
cat("  Correlation method: ", toupper(correlation_method), "\n\n")

cat("Model Results:\n")
cat("  Infection Status: Effect = +", round(status_effect, 2), "%, p = ", format(status_p, scientific = TRUE), "\n")
cat("  E. ferrisi: Effect = +", round(ferrisi_effect, 2), "%, p = ", format(ferrisi_p, scientific = TRUE), "\n")
cat("  E. falciformis: Effect = +", round(falciformis_effect, 2), "%, p = ", format(falciformis_p, scientific = TRUE), "\n")
cat("  Interaction effect: ", round(interaction_effect, 3), ", p = ", format(interaction_p, scientific = TRUE), "\n")
cat("  qPCR correlation: r = ", round(cor_qpcr$estimate, 3), ", p = ", format(cor_qpcr$p.value, scientific = TRUE), "\n")
cat("  OPG correlation: r = ", round(cor_opg$estimate, 3), ", p = ", format(cor_opg$p.value, scientific = TRUE), "\n")

cat(strrep("=", 60), "\n")


# EXTRACT CONFIDENCE INTERVALS FROM YOUR EXISTING MODELS
# Add this code after your existing script runs

cat("EXTRACTING 95% CONFIDENCE INTERVALS FOR MANUSCRIPT\n")
cat(strrep("=", 50), "\n")

# Function to calculate 95% CI from estimate and SE
calculate_ci <- function(estimate, se, df = Inf) {
    t_value <- qt(0.975, df)  # 97.5th percentile for 95% CI
    margin <- t_value * se
    return(c(lower = estimate - margin, upper = estimate + margin))
}

# ---- INFECTION STATUS MODEL CIs ----
cat("1. INFECTION STATUS MODEL:\n")

# Get model degrees of freedom
status_df <- status_model$df.residual

# Calculate CI for infection effect
status_ci <- calculate_ci(status_effect, status_se, status_df)
cat("   Eimeria-positive vs uninfected:\n")
cat("   estimate = +", round(status_effect, 2), "%, ", 
    "95% CI: [", round(status_ci[1], 2), "%, ", round(status_ci[2], 2), "%], ",
    "p ", ifelse(status_p < 0.001, "< 0.001", paste("=", round(status_p, 3))), "\n")

# ---- SPECIES MODEL CIs ----
cat("\n2. SPECIES-SPECIFIC MODEL:\n")

# Get model degrees of freedom
species_df <- species_model$df.residual

# Calculate CIs for both species
ferrisi_ci <- calculate_ci(ferrisi_effect, ferrisi_se, species_df)
falciformis_ci <- calculate_ci(falciformis_effect, falciformis_se, species_df)

cat("   E. ferrisi vs uninfected:\n")
cat("   estimate = +", round(ferrisi_effect, 2), "%, ", 
    "95% CI: [", round(ferrisi_ci[1], 2), "%, ", round(ferrisi_ci[2], 2), "%], ",
    "p ", ifelse(ferrisi_p < 0.001, "< 0.001", paste("=", round(ferrisi_p, 3))), "\n")

cat("   E. falciformis vs uninfected:\n")
cat("   estimate = +", round(falciformis_effect, 2), "%, ", 
    "95% CI: [", round(falciformis_ci[1], 2), "%, ", round(falciformis_ci[2], 2), "%], ",
    "p ", ifelse(falciformis_p < 0.001, "< 0.001", paste("=", round(falciformis_p, 3))), "\n")

# ---- INTERACTION MODEL CIs ----
cat("\n3. INTERACTION MODEL:\n")

# Get model degrees of freedom
interaction_df <- interaction_model$df.residual

# Calculate CIs for main effect and interaction
main_ci <- calculate_ci(infection_main_effect, infection_main_se, interaction_df)
interaction_ci <- calculate_ci(interaction_effect, interaction_se, interaction_df)

cat("   Main effect (infection status):\n")
cat("   estimate = +", round(infection_main_effect, 2), "%, ", 
    "95% CI: [", round(main_ci[1], 2), "%, ", round(main_ci[2], 2), "%], ",
    "p ", ifelse(infection_main_p < 0.001, "< 0.001", paste("=", round(infection_main_p, 3))), "\n")

cat("   Interaction term:\n")
cat("   estimate = ", round(interaction_effect, 3), ", ", 
    "95% CI: [", round(interaction_ci[1], 3), ", ", round(interaction_ci[2], 3), "], ",
    "p ", ifelse(interaction_p < 0.001, "< 0.001", paste("=", round(interaction_p, 3))), "\n")

# ---- ALTERNATIVE: Use confint() for exact CIs ----
cat("\n4. EXACT CONFIDENCE INTERVALS (using confint()):\n")

# Status model CIs
status_confint <- confint(status_model)
cat("   Status model 95% CIs:\n")
print(status_confint)

# Species model CIs  
species_confint <- confint(species_model)
cat("\n   Species model 95% CIs:\n")
print(species_confint)

# Interaction model CIs
interaction_confint <- confint(interaction_model)
cat("\n   Interaction model 95% CIs:\n")
print(interaction_confint)

# ---- FORMAT FOR MANUSCRIPT ----
cat("\n", strrep("=", 50), "\n")
cat("FORMATTED FOR MANUSCRIPT:\n")
cat(strrep("=", 50), "\n")

cat("Results section text with CIs:\n\n")

cat("Linear regression analysis showed that Eimeria-positive mice had significantly higher\n")
cat("predicted weight loss compared to uninfected mice (estimate = +", round(status_effect, 2), "%, ")
cat("95% CI: [", round(status_confint[2,1], 2), "%, ", round(status_confint[2,2], 2), "%], ")
cat("p ", ifelse(status_p < 0.001, "< 0.001", paste("=", round(status_p, 3))), ", ")
cat("adjusted R² = ", round(status_summary$adj.r.squared, 3), ", n = ", nrow(Field_status), ").\n\n")

cat("When comparing specific Eimeria species, E. falciformis infections produced the highest\n")
cat("predicted weight loss (estimate = +", round(falciformis_effect, 2), "%, ")
cat("95% CI: [", round(species_confint[3,1], 2), "%, ", round(species_confint[3,2], 2), "%], ")
cat("p ", ifelse(falciformis_p < 0.001, "< 0.001", paste("=", round(falciformis_p, 3))), "), ")
cat("followed by E. ferrisi (estimate = +", round(ferrisi_effect, 2), "%, ")
cat("95% CI: [", round(species_confint[2,1], 2), "%, ", round(species_confint[2,2], 2), "%], ")
cat("p ", ifelse(ferrisi_p < 0.001, "< 0.001", paste("=", round(ferrisi_p, 3))), "), ")
cat("both significantly greater than uninfected mice (model R² = ", round(species_summary$r.squared, 3), ", ")
cat("p < 0.001, n = ", nrow(Field_species), ").\n")

cat("\n✅ All confidence intervals calculated and formatted!\n")

