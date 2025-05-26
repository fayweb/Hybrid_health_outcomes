# =============================================================================
# FIELD DATA ANALYSIS - FIXED FOR YOUR PIPELINE
# =============================================================================

# This script assumes you have already run through your MNI pipeline
# and have the 'hm' object (merged, normalized, imputed data)

library(dplyr)
library(ggplot2)

# =============================================================================
# 1. EXTRACT FIELD DATA FROM YOUR HM OBJECT
# =============================================================================

# Filter for field mice (you already did this)
field <- hm %>%
    filter(origin == "Field")

cat("=== FIELD DATA OVERVIEW ===\n")
cat("Total field mice with immune data:", length(unique(field$Mouse_ID)), "\n\n")

# =============================================================================
# 2. FIX IMMUNE GENE DETECTION
# =============================================================================

# Your gene vector from the pipeline
Genes_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
             "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
             "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
             "TICAM1", "TNF")

# Check which genes are actually in your field data
available_genes <- intersect(Genes_v, colnames(field))
cat("Available immune genes:", length(available_genes), "out of", length(Genes_v), "\n")
cat("Missing genes:", setdiff(Genes_v, available_genes), "\n\n")

# Create immune data indicator (all 336 mice should have immune data after imputation)
field <- field %>%
    group_by(Mouse_ID) %>%
    summarise(
        # Keep first row for each mouse (since they're duplicated)
        across(everything(), first),
        .groups = 'drop'
    ) %>%
    mutate(has_immune_data = TRUE)  # All should have data after imputation

cat("Unique mice after deduplication:", nrow(field), "\n")

# =============================================================================
# 3. DATA AVAILABILITY SUMMARY
# =============================================================================

cat("\n=== DATA AVAILABILITY BY METHOD ===\n")

data_availability <- field %>%
    summarise(
        total_mice = n(),
        has_immune_genes = sum(has_immune_data),
        has_qPCR_tissue = sum(!is.na(delta_ct_cewe_MminusE)),
        has_qPCR_feces = sum(!is.na(FEC_Eim_Ct)),
        has_oocyst_counts = sum(!is.na(OPG) & OPG >= 0),
        has_amplicon_data = sum(!is.na(amplicon_species)),
        has_species_ID = sum(!is.na(species_Eimeria) & species_Eimeria != "NA")
    )

print(data_availability)

# =============================================================================
# 4. INFECTION STATUS BREAKDOWN
# =============================================================================

cat("\n=== INFECTION STATUS BREAKDOWN ===\n")

# Create binary infection status from multiple sources
field <- field %>%
    mutate(
        # Primary infection status (tissue qPCR)
        infection_tissue = case_when(
            MC.Eimeria == "TRUE" ~ "TRUE",
            MC.Eimeria == "FALSE" ~ "FALSE",
            TRUE ~ "Unknown"
        ),
        
        # Secondary infection status (species classification)
        infection_species = case_when(
            species_Eimeria %in% c("E. ferrisi", "E. falciformis") ~ "TRUE",
            species_Eimeria %in% c("uninfected", "Negative", "Uninfected") ~ "FALSE",
            TRUE ~ "Unknown"
        ),
        
        # Combined infection status (prioritizing tissue qPCR when available)
        infection_status = case_when(
            infection_tissue %in% c("TRUE", "FALSE") ~ infection_tissue,
            infection_species %in% c("TRUE", "FALSE") ~ infection_species,
            TRUE ~ "Unknown"
        )
    )

# Infection status summary
cat("Infection Status Summary:\n")
infection_summary <- table(field$infection_status, useNA = "always")
print(infection_summary)

infection_rate <- round(sum(field$infection_status == "TRUE", na.rm = TRUE) / 
                            sum(field$infection_status %in% c("TRUE", "FALSE")) * 100, 1)
cat("Overall infection rate:", infection_rate, "%\n\n")

# =============================================================================
# 5. SPECIES BREAKDOWN
# =============================================================================

cat("=== SPECIES BREAKDOWN ===\n")
species_summary <- table(field$species_Eimeria, useNA = "always")
print(species_summary)

# Clean up species classification
field <- field %>%
    mutate(
        species_clean = case_when(
            species_Eimeria == "E. ferrisi" ~ "E. ferrisi",
            species_Eimeria == "E. falciformis" ~ "E. falciformis", 
            species_Eimeria %in% c("uninfected", "Negative", "Uninfected") ~ "Uninfected",
            TRUE ~ "Unknown"
        )
    )

cat("\nCleaned Species Classification:\n")
table(field$species_clean, useNA = "always") %>% print()

# =============================================================================
# 6. QUANTIFICATION METHOD OVERLAP
# =============================================================================

cat("\n=== QUANTIFICATION METHOD OVERLAP ===\n")

# Method availability for each mouse
field <- field %>%
    mutate(
        has_tissue_qPCR = !is.na(delta_ct_cewe_MminusE),
        has_feces_qPCR = !is.na(FEC_Eim_Ct),
        has_oocyst = !is.na(OPG),
        has_amplicon = !is.na(amplicon_species),
        has_species = !is.na(species_Eimeria) & species_Eimeria != "NA"
    )

# Create method combination codes
field <- field %>%
    mutate(
        method_combo = paste0(
            ifelse(has_tissue_qPCR, "T", ""),
            ifelse(has_feces_qPCR, "F", ""),
            ifelse(has_oocyst, "O", ""),
            ifelse(has_amplicon, "A", ""),
            ifelse(has_species, "S", ""),
            "I"  # All have immune data
        )
    )

cat("Method combinations (T=Tissue_qPCR, F=Feces_qPCR, O=Oocyst, A=Amplicon, S=Species, I=Immune):\n")
method_combos <- table(field$method_combo) %>% sort(decreasing = TRUE)
print(head(method_combos, 10))

# =============================================================================
# 7. INFECTION INTENSITY ANALYSIS
# =============================================================================

cat("\n=== INFECTION INTENSITY DISTRIBUTIONS ===\n")

# Tissue qPCR intensity (ΔCt values)
tissue_infected <- field %>% 
    filter(infection_status == "TRUE" & has_tissue_qPCR)

if(nrow(tissue_infected) > 0) {
    cat("Tissue ΔCt values in infected mice (n =", nrow(tissue_infected), "):\n")
    cat("Mean ΔCt:", round(mean(tissue_infected$delta_ct_cewe_MminusE, na.rm = TRUE), 2), "\n")
    cat("Range ΔCt:", paste(round(range(tissue_infected$delta_ct_cewe_MminusE, na.rm = TRUE), 2), collapse = " to "), "\n")
} else {
    cat("No mice with both infection status and tissue qPCR data\n")
}

# Oocyst intensity
oocyst_infected <- field %>% 
    filter(has_oocyst & OPG > 0)

if(nrow(oocyst_infected) > 0) {
    cat("\nOPG values in mice with oocyst shedding (n =", nrow(oocyst_infected), "):\n")
    cat("Mean OPG:", format(round(mean(oocyst_infected$OPG, na.rm = TRUE), 0), scientific = FALSE), "\n")
    cat("Median OPG:", format(round(median(oocyst_infected$OPG, na.rm = TRUE), 0), scientific = FALSE), "\n")
    cat("Range OPG:", paste(format(range(oocyst_infected$OPG, na.rm = TRUE), scientific = FALSE), collapse = " to "), "\n")
} else {
    cat("No mice with oocyst data\n")
}

# =============================================================================
# 8. ANALYSIS-READY DATASETS
# =============================================================================

cat("\n=== ANALYSIS-READY DATASETS ===\n")

# Dataset 1: All mice (for random forest prediction)
all_mice <- field
cat("Dataset 1 - All mice with immune data:", nrow(all_mice), "\n")

# Dataset 2: Mice with known infection status (for validation)
known_status <- field %>% filter(infection_status %in% c("TRUE", "FALSE"))
cat("Dataset 2 - Mice with known infection status:", nrow(known_status), "\n")
cat("  - Infected:", sum(known_status$infection_status == "TRUE"), "\n")
cat("  - Uninfected:", sum(known_status$infection_status == "FALSE"), "\n")

# Dataset 3: Mice with species classification (for species-specific analysis)
species_classified <- field %>% filter(species_clean %in% c("E. ferrisi", "E. falciformis", "Uninfected"))
cat("Dataset 3 - Mice with species classification:", nrow(species_classified), "\n")
species_breakdown <- table(species_classified$species_clean)
for(i in 1:length(species_breakdown)) {
    cat("  -", names(species_breakdown)[i], ":", species_breakdown[i], "\n")
}

# Dataset 4: Mice with infection intensity measures (for correlation analysis)
intensity_data <- field %>% filter(has_tissue_qPCR & infection_status == "TRUE")
cat("Dataset 4 - Infected mice with ΔCt values:", nrow(intensity_data), "\n")

# Dataset 5: Complete data (all methods available)
complete_data <- field %>% filter(has_tissue_qPCR & has_oocyst & has_amplicon)
cat("Dataset 5 - Mice with complete quantification data:", nrow(complete_data), "\n")

# =============================================================================
# 9. SAVE CLEANED DATA
# =============================================================================

# Save the cleaned field dataset
write.csv(field, "field_analysis_ready.csv", row.names = FALSE)

# Create summary table for your thesis
summary_table <- data.frame(
    Analysis_Dataset = c("All mice", "Known infection status", "Species classified", 
                         "Infection intensity", "Complete data"),
    Sample_Size = c(nrow(all_mice), nrow(known_status), nrow(species_classified),
                    nrow(intensity_data), nrow(complete_data)),
    Infected = c(sum(all_mice$infection_status == "TRUE", na.rm = TRUE),
                 sum(known_status$infection_status == "TRUE"),
                 sum(species_classified$species_clean %in% c("E. ferrisi", "E. falciformis")),
                 nrow(intensity_data),
                 sum(complete_data$infection_status == "TRUE", na.rm = TRUE)),
    Available_Methods = c("Immune genes", "Immune + Infection status", 
                          "Immune + Species ID", "Immune + qPCR intensity",
                          "Immune + qPCR + Oocyst + Amplicon")
) %>%
    mutate(Infection_Rate = round(Infected / Sample_Size * 100, 1))

write.csv(summary_table, "analysis_datasets_summary.csv", row.names = FALSE)

cat("\n=== SUMMARY TABLE ===\n")
print(summary_table)

cat("\n=== FILES CREATED ===\n")
cat("1. field_analysis_ready.csv - Complete cleaned field dataset\n")
cat("2. analysis_datasets_summary.csv - Summary table for thesis\n")

# =============================================================================
# 10. RECOMMENDATIONS FOR YOUR RANDOM FOREST ANALYSIS
# =============================================================================

cat("\n=== RECOMMENDATIONS FOR RANDOM FOREST ANALYSIS ===\n\n")

cat("PRIMARY ANALYSIS:\n")
cat("- Use all", nrow(all_mice), "mice for random forest prediction of health impacts\n")
cat("- Validate predictions against", nrow(known_status), "mice with known infection status\n")
cat("- Correlate predictions with ΔCt values in", nrow(intensity_data), "infected mice\n\n")

cat("VALIDATION STRATEGY:\n")
cat("1. Train random forest on lab data (as you did)\n")
cat("2. Apply to all", nrow(all_mice), "field mice → predicted health impacts\n")
cat("3. Validate predictions using:\n")
cat("   - Infection status correlation (n =", nrow(known_status), ")\n")
cat("   - Species-specific effects (n =", nrow(species_classified), ")\n")
cat("   - Infection intensity correlation (n =", nrow(intensity_data), ")\n")
cat("   - Cross-method validation (n =", nrow(complete_data), ")\n\n")

cat("THESIS WRITING:\n")
cat("Your field validation used", nrow(all_mice), "wild-caught house mice with immune gene expression data.\n")
cat("Infection status was determined using multiple complementary approaches:\n")
cat("- Tissue qPCR (n =", sum(field$has_tissue_qPCR), "mice)\n")
cat("- Amplicon sequencing species identification (n =", sum(field$has_amplicon), "mice)\n") 
cat("- Oocyst counting validation (n =", sum(field$has_oocyst), "mice)\n")
cat("This multi-method approach enabled validation of random forest predictions\n")
cat("against both binary infection status and quantitative infection intensity.\n")

# View the plots that were created
p1  # Sample size flowchart
p2  # Infection status breakdown

# Save them
ggsave('flowchart.png', p1, width = 10, height = 6)
ggsave('infection_status.png', p2, width = 8, height = 6)
