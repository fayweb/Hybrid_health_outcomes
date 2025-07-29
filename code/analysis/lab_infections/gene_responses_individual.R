# ***********************************************************
# SCRIPT 1: Cross-Population Individual Gene Analysis
# Purpose: Generate Figure 3A-B and specific estimates reported in results
# Results: Lab vs Field gene responses, identify cross-validated genes
# ***********************************************************

library(tidyverse)
library(broom)
library(ggplot2)
library(patchwork)

# Define immune gene panel (19 genes)
Genes_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", 
             "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", 
             "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

# ***********************************************************
# DATA PREPARATION
# ***********************************************************

# Laboratory mice (n=136)
lab_data <- hm %>% 
    filter(origin == "Lab") %>%
    dplyr::select(Mouse_ID, current_infection, all_of(Genes_v))

# Field mice with species data (n=169)
field_data <- hm %>% 
    filter(origin == "Field") %>%
    drop_na(species_Eimeria) %>%
    dplyr::select(Mouse_ID, species_Eimeria, all_of(Genes_v))

cat("Sample sizes: Lab n =", nrow(lab_data), ", Field n =", nrow(field_data), "\n")

# ***********************************************************
# LABORATORY REGRESSION ANALYSES
# ***********************************************************

# Run regression for each gene in laboratory
lab_results <- map(Genes_v, function(gene) {
    formula_str <- paste(gene, "~ current_infection")
    lm(as.formula(formula_str), data = lab_data)
}) %>% set_names(Genes_v)

# Extract coefficients with confidence intervals
lab_coefficients <- map_dfr(names(lab_results), function(gene) {
    tidy(lab_results[[gene]], conf.int = TRUE) %>%
        mutate(gene = gene, population = "Laboratory")
}) %>%
    filter(!str_detect(term, "Intercept")) %>%
    mutate(
        species = case_when(
            str_detect(term, "E. falciformis") ~ "E. falciformis",
            str_detect(term, "E. ferrisi") ~ "E. ferrisi",
            TRUE ~ term
        )
    ) %>%
    dplyr::select(gene, species, estimate, conf.low, conf.high, p.value, population)

# ***********************************************************
# FIELD REGRESSION ANALYSES  
# ***********************************************************

# Run regression for each gene in field
field_results <- map(Genes_v, function(gene) {
    formula_str <- paste(gene, "~ species_Eimeria")
    lm(as.formula(formula_str), data = field_data)
}) %>% set_names(Genes_v)

# Extract coefficients
field_coefficients <- map_dfr(names(field_results), function(gene) {
    tidy(field_results[[gene]], conf.int = TRUE) %>%
        mutate(gene = gene, population = "Field")
}) %>%
    filter(!str_detect(term, "Intercept")) %>%
    mutate(
        species = case_when(
            str_detect(term, "E. falciformis") ~ "E. falciformis", 
            str_detect(term, "E. ferrisi") ~ "E. ferrisi",
            TRUE ~ term
        )
    ) %>%
    dplyr::select(gene, species, estimate, conf.low, conf.high, p.value, population)

# ***********************************************************
# IDENTIFY CROSS-VALIDATED GENES
# ***********************************************************

# Combine results and find genes significant in both populations
combined_results <- bind_rows(lab_coefficients, field_coefficients)

cross_validated <- combined_results %>%
    filter(p.value < 0.05) %>%
    count(gene, species) %>%
    filter(n == 2) %>%  # Significant in both populations
    arrange(gene, species)

cat("\n=== CROSS-POPULATION VALIDATED GENES ===\n")
print(cross_validated)

# Extract specific estimates mentioned in results
specific_estimates <- combined_results %>%
    filter(gene %in% c("CXCL9", "IFNy", "IDO1", "TICAM1", "NCR1", "PRF1")) %>%
    arrange(gene, species, population)

cat("\n=== KEY ESTIMATES FOR MANUSCRIPT ===\n")
print(specific_estimates, n = Inf)

# ***********************************************************
# CREATE FIGURE 3A-B
# ***********************************************************

# Color scheme
color_scheme <- c("E. falciformis" = "#FF0000", "E. ferrisi" = "#7A0092")

# Panel A: Laboratory results
panel_A <- ggplot(lab_coefficients, aes(x = gene, y = estimate, color = species)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0.2, position = position_dodge(0.5)) +
    geom_point(size = 2, position = position_dodge(0.5)) +
    scale_color_manual(values = color_scheme, name = "Species") +
    labs(title = "A) Laboratory infections",
         x = "Immune genes", 
         y = "Regression coefficient (95% CI)") +
    theme_minimal() +
    theme(legend.position = "top")

# Panel B: Field results
panel_B <- ggplot(field_coefficients, aes(x = gene, y = estimate, color = species)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0.2, position = position_dodge(0.5)) +
    geom_point(size = 2, position = position_dodge(0.5)) +
    scale_color_manual(values = color_scheme, name = "Species") +
    labs(title = "B) Field infections",
         x = "Immune genes", 
         y = "Regression coefficient (95% CI)") +
    theme_minimal() +
    theme(legend.position = "none")

# Combine panels
figure_3AB <- panel_A | panel_B

print(figure_3AB)

# Save figure
ggsave(paste0(an_fi, "/figure_3AB_gene_responses.pdf"), 
       figure_3AB, width = 12, height = 8, dpi = 300)

# ***********************************************************
# EXPORT RESULTS FOR SUPPLEMENTARY TABLE 5
# ***********************************************************

supplementary_table5 <- combined_results %>%
    mutate(
        ci_text = paste0("(", round(conf.low, 2), " to ", round(conf.high, 2), ")"),
        estimate_text = paste0(round(estimate, 2), " ", ci_text),
        significance = ifelse(p.value < 0.05, "Yes", "No")
    ) %>%
    dplyr::select(gene, species, population, estimate_text, p.value, significance) %>%
    pivot_wider(names_from = population, 
                values_from = c(estimate_text, p.value, significance),
                names_sep = "_")

write_csv(supplementary_table5, paste0(tables, "/supplementary_table5.csv"))

cat("\n✅ Analysis complete! Cross-validated genes:", nrow(cross_validated), "\n")
cat("📊 Figure 3A-B saved, Supplementary Table 5 exported\n")


#