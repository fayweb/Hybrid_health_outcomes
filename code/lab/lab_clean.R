# ============================================================================
# SCRIPT: lab_clean.R
# PURPOSE: Clean and standardize laboratory infection data
# AUTHOR: Fay Webster
# 
# INPUTS: 
#   - lab (dataframe from lab_import.R)
#   - Path variables (dlab_final)
#
# OUTPUTS:
#   - Challenge (cleaned lab dataframe)
#   - lab_cleaned_data.csv (saved to dlab_final)
#
# NOTES: This script performs data standardization and creates infection
#        categories. Visualization is handled separately in lab_visualize.R
# ============================================================================

cat("Starting laboratory data cleaning...\n")

# Validate inputs
if (!exists("lab")) {
    stop("Error: 'lab' dataframe not found. Please run lab_import.R first.")
}

# ============================================================================
# STEP 1: STANDARDIZE PARASITE STRAIN NAMES
# ============================================================================

cat("1. Standardizing parasite strain names...\n")

# Create standardized parasite names for primary infections
lab <- lab %>%
    dplyr::mutate(
        Parasite_primary = case_when(
            primary_infection == "E64" ~ "E_ferrisi",
            primary_infection == "E88" ~ "E_falciformis", 
            primary_infection == "Eflab" ~ "E_falciformis",
            primary_infection == "E139" ~ "E_ferrisi",
            primary_infection == "UNI" ~ "uninfected",
            TRUE ~ NA_character_  # Use NA instead of empty string
        )
    )

# Create standardized parasite names for challenge infections  
lab <- lab %>%
    dplyr::mutate(
        Parasite_challenge = case_when(
            challenge_infection == "E64" ~ "E_ferrisi",
            challenge_infection == "E88" ~ "E_falciformis",
            challenge_infection == "Eflab" ~ "E_falciformis", 
            challenge_infection == "E139" ~ "E_ferrisi",
            challenge_infection == "UNI" ~ "uninfected",
            TRUE ~ NA_character_
        )
    )

# Check for any unmapped infections
unmapped_primary <- lab %>% 
    filter(is.na(Parasite_primary)) %>% 
    select(primary_infection) %>% 
    distinct()

unmapped_challenge <- lab %>%
    filter(is.na(Parasite_challenge)) %>%
    select(challenge_infection) %>%
    distinct()

if (nrow(unmapped_primary) > 0) {
    warning("Unmapped primary infections found: ", 
            paste(unmapped_primary$primary_infection, collapse = ", "))
}

if (nrow(unmapped_challenge) > 0) {
    warning("Unmapped challenge infections found: ",
            paste(unmapped_challenge$challenge_infection, collapse = ", "))
}

# ============================================================================
# STEP 2: CREATE INFECTION HISTORY CATEGORIES
# ============================================================================

cat("2. Creating infection history categories...\n")

lab <- lab %>%
    dplyr::mutate(
        infection_history = case_when(
            Parasite_primary == "uninfected" & Parasite_challenge == "uninfected" ~ "uninfected",
            Parasite_primary == "uninfected" & Parasite_challenge == "E_ferrisi" ~ "uninfected_ferrisi", 
            Parasite_primary == "uninfected" & Parasite_challenge == "E_falciformis" ~ "uninfected_falciformis",
            Parasite_primary == "E_falciformis" & Parasite_challenge == "E_falciformis" ~ "falciformis_falciformis",
            Parasite_primary == "E_falciformis" & Parasite_challenge == "E_ferrisi" ~ "falciformis_ferrisi",
            Parasite_primary == "E_falciformis" & Parasite_challenge == "uninfected" ~ "falciformis_uninfected",
            Parasite_primary == "E_ferrisi" & Parasite_challenge == "E_falciformis" ~ "ferrisi_falciformis",
            Parasite_primary == "E_ferrisi" & Parasite_challenge == "E_ferrisi" ~ "ferrisi_ferrisi",
            Parasite_primary == "E_ferrisi" & Parasite_challenge == "uninfected" ~ "ferrisi_uninfected",
            TRUE ~ NA_character_
        )
    )

# Check for any unmapped infection histories
unmapped_history <- lab %>%
    filter(is.na(infection_history)) %>%
    select(Parasite_primary, Parasite_challenge) %>%
    distinct()

if (nrow(unmapped_history) > 0) {
    warning("Unmapped infection histories found:")
    print(unmapped_history)
}

# ============================================================================
# STEP 3: DATA QUALITY CONTROL
# ============================================================================

cat("3. Performing data quality control...\n")

# Replace infinite values with NA
lab[sapply(lab, is.infinite)] <- NA

# Clean weight data
# Replace 0 weights with NA (incorrect measurements)
weight_zeros <- sum(lab$weight == 0, na.rm = TRUE)
if (weight_zeros > 0) {
    cat("  - Replacing", weight_zeros, "zero weight values with NA\n")
    lab$weight[lab$weight == 0] <- NA
}

# Replace 0 relative weights with NA
rel_weight_zeros <- sum(lab$relative_weight == 0, na.rm = TRUE)
if (rel_weight_zeros > 0) {
    cat("  - Replacing", rel_weight_zeros, "zero relative weight values with NA\n")
    lab$relative_weight[lab$relative_weight == 0] <- NA
}

# Set dpi to NA when weight is missing (after death)
dpi_na_count <- sum(is.na(lab$weight))
lab$dpi[is.na(lab$weight)] <- NA
cat("  - Set", dpi_na_count, "dpi values to NA (post-death timepoints)\n")

# ============================================================================
# STEP 4: CALCULATE DERIVED VARIABLES
# ============================================================================

cat("4. Calculating derived variables...\n")

# Calculate maximum DPI and weight loss for each mouse
lab <- lab %>%
    dplyr::filter(!is.na(weight)) %>%  # Remove NA weights for calculation
    dplyr::group_by(EH_ID, infection) %>%
    dplyr::mutate(
        max_dpi = max(dpi, na.rm = TRUE),
        WL_max = 100 - min(relative_weight, na.rm = TRUE)  # Maximum weight loss
    ) %>%
    dplyr::ungroup()

cat("  - Calculated max_dpi and WL_max for", 
    length(unique(lab$EH_ID)), "mice\n")

# ============================================================================
# STEP 5: GENE EXPRESSION CLEANUP
# ============================================================================

cat("5. Cleaning gene expression data...\n")

# Handle CXCR3 columns (keep the bio version)
if ("CXCR3" %in% colnames(lab) && "CXCR3_bio" %in% colnames(lab)) {
    cat("  - Removing duplicate CXCR3 column, keeping CXCR3_bio\n")
    lab <- lab %>% 
        dplyr::select(-CXCR3) %>%
        dplyr::rename(CXCR3 = CXCR3_bio)
} else if ("CXCR3_bio" %in% colnames(lab)) {
    cat("  - Renaming CXCR3_bio to CXCR3\n")
    lab <- lab %>%
        dplyr::rename(CXCR3 = CXCR3_bio)
}

# ============================================================================
# STEP 6: STANDARDIZE COLUMN NAMES
# ============================================================================

cat("6. Standardizing column names...\n")

# Add origin identifier
lab <- lab %>%
    dplyr::mutate(origin = "Lab")

# Rename columns to match field data
lab <- lab %>% 
    dplyr::rename(
        Mouse_ID = EH_ID,
        delta_ct_cewe_MminusE = delta,
        MC.Eimeria = Eim_MC,
        Feces_Weight = feces_weight
    )

cat("  - Renamed key columns for consistency\n")

# ============================================================================
# STEP 7: DETERMINE CURRENT INFECTION STATUS
# ============================================================================

cat("7. Determining current infection status...\n")

# Create current infection status based on challenge infection and qPCR results
lab <- lab %>%
    dplyr::mutate(
        current_infection = case_when(
            Parasite_challenge == "E_ferrisi" & MC.Eimeria == TRUE ~ "E_ferrisi",
            Parasite_challenge == "E_ferrisi" & MC.Eimeria == FALSE ~ "uninfected",
            Parasite_challenge == "E_falciformis" & MC.Eimeria == TRUE ~ "E_falciformis", 
            Parasite_challenge == "E_falciformis" & MC.Eimeria == FALSE ~ "uninfected",
            Parasite_challenge == "uninfected" & MC.Eimeria == TRUE ~ "E_falciformis",  # Unexpected infection
            Parasite_challenge == "uninfected" & MC.Eimeria == FALSE ~ "uninfected",
            TRUE ~ Parasite_challenge
        )
    )

# ============================================================================
# STEP 8: CREATE IMMUNIZATION CATEGORIES
# ============================================================================

cat("8. Creating immunization categories...\n")

lab <- lab %>%
    dplyr::mutate(
        immunization = case_when(
            infection_history == "falciformis_ferrisi" ~ "heterologous",
            infection_history == "ferrisi_falciformis" ~ "heterologous", 
            infection_history == "falciformis_uninfected" ~ "primary_only",
            infection_history == "ferrisi_uninfected" ~ "primary_only",
            infection_history == "ferrisi_ferrisi" ~ "homologous",
            infection_history == "falciformis_falciformis" ~ "homologous",
            infection_history == "uninfected_falciformis" ~ "naive",
            infection_history == "uninfected_ferrisi" ~ "naive", 
            infection_history == "uninfected" ~ "uninfected",
            TRUE ~ NA_character_
        )
    )

# ============================================================================
# STEP 9: DATA VALIDATION AND SUMMARY
# ============================================================================

cat("9. Validating cleaned data...\n")

# Check sample sizes by group
sample_summary <- lab %>%
    group_by(current_infection, immunization) %>%
    summarise(
        n_mice = n_distinct(Mouse_ID),
        .groups = "drop"
    )

cat("  Sample sizes by group:\n")
print(sample_summary)

# Check for missing critical variables
critical_vars <- c("Mouse_ID", "current_infection", "WL_max", "MC.Eimeria")
missing_summary <- lab %>%
    summarise(across(all_of(critical_vars), ~sum(is.na(.))))

cat("  Missing values in critical variables:\n")
print(missing_summary)

# ============================================================================
# STEP 10: SAVE CLEANED DATA
# ============================================================================

cat("10. Saving cleaned data...\n")

# Create final cleaned dataset
Challenge <- lab

# Ensure output directory exists
if (!dir.exists(dlab_final)) {
    dir.create(dlab_final, recursive = TRUE)
    cat("  - Created output directory:", dlab_final, "\n")
}

# Save cleaned data
output_file <- file.path(dlab_final, "lab_cleaned_data.csv")
write.csv(Challenge, output_file, row.names = FALSE)

cat("  - Cleaned data saved to:", output_file, "\n")

# ============================================================================
# COMPLETION SUMMARY
# ============================================================================

cat("\n=== LAB DATA CLEANING COMPLETED ===\n")
cat("Final dataset summary:\n")
cat("  - Total observations:", nrow(Challenge), "\n")
cat("  - Unique mice:", length(unique(Challenge$Mouse_ID)), "\n")
cat("  - Infection groups:", length(unique(Challenge$current_infection)), "\n")
cat("  - Immunization categories:", length(unique(Challenge$immunization)), "\n")

# Display infection group summary
infection_summary <- Challenge %>%
    count(current_infection, name = "n_observations") %>%
    arrange(desc(n_observations))

cat("  - Observations by infection group:\n")
print(infection_summary)

cat("\nReady for next step: lab_visualize.R\n")
cat("======================================\n")
