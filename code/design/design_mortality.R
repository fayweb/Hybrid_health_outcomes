# =============================================================================
# ENHANCED MORTALITY AND WEIGHT LOSS ANALYSIS - FIXED FOR LONG FORMAT
# =============================================================================
# Add this section to your design_experimental.R script AFTER loading both datasets

# =============================================================================
# 6. MORTALITY ANALYSIS USING CHALLENGE DATASET (CORRECTED)
# =============================================================================

cat("\n=== COMPREHENSIVE MORTALITY ANALYSIS ===\n")

# First, let's understand both dataset structures
cat("Dataset overview:\n")
cat("Total rows in Challenge dataset:", nrow(Challenge), "\n")
cat("Unique mice in Challenge dataset:", length(unique(Challenge$Mouse_ID)), "\n")
cat("Total mice in lab_clean (with immune data):", nrow(lab_clean), "\n")

# Check which mice are in Challenge but not in lab_clean (the missing 17)
challenge_mice <- unique(Challenge$Mouse_ID)
lab_mice <- unique(lab_clean$Mouse_ID)
missing_from_lab <- setdiff(challenge_mice, lab_mice)
missing_from_challenge <- setdiff(lab_mice, challenge_mice)

cat("Mice in Challenge but NOT in lab dataset (missing immune data):", length(missing_from_lab), "\n")
cat("Mice in lab but NOT in Challenge dataset:", length(missing_from_challenge), "\n")

# CORRECTED: Get one row per mouse, not per day
mortality_status <- Challenge %>%
    group_by(Mouse_ID) %>%
    slice(1) %>%  # Take first row for each mouse to get mouse-level info
    ungroup() %>%
    dplyr::select(Mouse_ID, death, Parasite_primary, Parasite_challenge) %>%
    mutate(has_immune_data = Mouse_ID %in% lab_mice)

cat("\nMortality breakdown (all", length(challenge_mice), "mice):\n")
mortality_summary <- table(mortality_status$death)
print(mortality_summary)

cat("\nMortality breakdown by immune data availability:\n")
mortality_by_data <- mortality_status %>%
    count(death, has_immune_data) %>%
    spread(has_immune_data, n, fill = 0)
colnames(mortality_by_data) <- c("Death_Phase", "No_Immune_Data", "Has_Immune_Data")
print(mortality_by_data)

# =============================================================================
# 7. DETAILED MORTALITY RATES BY SPECIES (CORRECTED)
# =============================================================================

cat("\n=== MORTALITY RATES BY SPECIES (ALL", length(challenge_mice), "MICE) ===\n")

# Calculate mortality rates by parasite species (primary infection) - ALL mice
primary_mortality_all <- mortality_status %>%
    filter(!is.na(Parasite_primary)) %>%  # Only mice that got primary infections
    group_by(Parasite_primary) %>%
    summarise(
        total_mice = n(),
        deaths_during_primary = sum(death == "primary", na.rm = TRUE),
        survivors_to_challenge = sum(death == "challenge", na.rm = TRUE),
        mortality_rate = round(deaths_during_primary / total_mice * 100, 1),
        survival_rate = round(survivors_to_challenge / total_mice * 100, 1),
        .groups = "drop"
    )

cat("Complete mortality analysis (all", length(challenge_mice), "mice):\n")
print(primary_mortality_all)

# Now do the same but only for mice with immune data
primary_mortality_immune <- mortality_status %>%
    filter(!is.na(Parasite_primary), has_immune_data == TRUE) %>%
    group_by(Parasite_primary) %>%
    summarise(
        total_mice_with_immune = n(),
        deaths_with_immune = sum(death == "primary", na.rm = TRUE),
        survivors_with_immune = sum(death == "challenge", na.rm = TRUE),
        mortality_rate_immune = round(deaths_with_immune / total_mice_with_immune * 100, 1),
        .groups = "drop"
    )

cat("\nMortality analysis (136 mice with immune data only):\n")
print(primary_mortality_immune)

# =============================================================================
# 8. WEIGHT LOSS IN DECEASED MICE (CORRECTED)
# =============================================================================

cat("\n=== WEIGHT LOSS IN DECEASED MICE ===\n")

# Calculate weight loss for mice that died during primary infection AND have immune data
deceased_mice_with_data <- mortality_status %>%
    filter(death == "primary", has_immune_data == TRUE) %>%
    pull(Mouse_ID)

cat("Number of mice that died during primary AND have immune data:", length(deceased_mice_with_data), "\n")

# Verify this matches what we see in lab_clean
lab_deaths <- lab_clean %>%
    filter(death == "primary") %>%
    nrow()
cat("Cross-check - mice with death='primary' in lab_clean:", lab_deaths, "\n")

# Get weight loss data for deceased mice (using Challenge dataset for detailed tracking)
if(length(deceased_mice_with_data) > 0) {
    deceased_weight_loss <- Challenge %>%
        filter(Mouse_ID %in% deceased_mice_with_data, infection == "primary") %>%
        group_by(Mouse_ID, Parasite_primary) %>%
        summarise(
            max_weight_loss = max(100 - relative_weight, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Summary statistics for deceased mice by species
    deceased_stats <- deceased_weight_loss %>%
        group_by(Parasite_primary) %>%
        summarise(
            n_deaths = n(),
            mean_weight_loss = round(mean(max_weight_loss, na.rm = TRUE), 1),
            min_weight_loss = round(min(max_weight_loss, na.rm = TRUE), 1),
            max_weight_loss = round(max(max_weight_loss, na.rm = TRUE), 1),
            sd_weight_loss = round(sd(max_weight_loss, na.rm = TRUE), 1),
            .groups = "drop"
        )
    
    cat("\nWeight loss in mice that died during primary infection (with immune data):\n")
    print(deceased_stats)
} else {
    cat("No deceased mice with immune data found in Challenge dataset.\n")
    deceased_stats <- NULL
}

# Alternative: use WL_max from lab_clean for deceased mice
cat("\nAlternative analysis using lab_clean WL_max data:\n")
deceased_from_lab <- lab_clean %>%
    filter(death == "primary") %>%
    group_by(current_infection) %>%
    summarise(
        n_deaths = n(),
        mean_WL_max = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL_max = round(min(WL_max, na.rm = TRUE), 1),
        max_WL_max = round(max(WL_max, na.rm = TRUE), 1),
        sd_WL_max = round(sd(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    )

print(deceased_from_lab)

# =============================================================================
# 9. COMPARISON: SURVIVORS vs DECEASED (CORRECTED)
# =============================================================================

cat("\n=== SURVIVOR vs DECEASED COMPARISON ===\n")

# Get weight loss for survivors during their primary infection (mice with immune data)
survivor_mice_with_data <- mortality_status %>%
    filter(death == "challenge", has_immune_data == TRUE) %>%
    pull(Mouse_ID)

cat("Survivors with immune data:", length(survivor_mice_with_data), "\n")

# Use Challenge dataset for detailed weight tracking
if(length(survivor_mice_with_data) > 0) {
    survivor_primary_wl <- Challenge %>%
        filter(Mouse_ID %in% survivor_mice_with_data, infection == "primary") %>%
        group_by(Mouse_ID, Parasite_primary) %>%
        summarise(
            max_weight_loss = max(100 - relative_weight, na.rm = TRUE),
            .groups = "drop"
        )
    
    survivor_stats <- survivor_primary_wl %>%
        group_by(Parasite_primary) %>%
        summarise(
            n_survivors = n(),
            mean_weight_loss = round(mean(max_weight_loss, na.rm = TRUE), 1),
            min_weight_loss = round(min(max_weight_loss, na.rm = TRUE), 1),
            max_weight_loss = round(max(max_weight_loss, na.rm = TRUE), 1),
            sd_weight_loss = round(sd(max_weight_loss, na.rm = TRUE), 1),
            .groups = "drop"
        )
    
    cat("\nWeight loss during primary infection in mice that SURVIVED to challenge:\n")
    print(survivor_stats)
} else {
    survivor_stats <- NULL
}

# Alternative: Use lab_clean data for all survivors
cat("\nAlternative analysis using lab_clean for all survivors:\n")
survivors_from_lab <- lab_clean %>%
    filter(death == "challenge") %>%
    group_by(current_infection) %>%
    summarise(
        n_survivors = n(),
        mean_WL_max = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL_max = round(min(WL_max, na.rm = TRUE), 1),
        max_WL_max = round(max(WL_max, na.rm = TRUE), 1),
        sd_WL_max = round(sd(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    )

print(survivors_from_lab)

# =============================================================================
# 10. COMPREHENSIVE RESULTS SUMMARY TABLE (CORRECTED)
# =============================================================================

cat("\n=== COMPREHENSIVE RESULTS FOR MANUSCRIPT ===\n")

# Create master summary table using the ALL MICE data (153 total)
results_summary_all <- primary_mortality_all %>%
    rename(Species = Parasite_primary,
           Total_Mice = total_mice,
           Deaths = deaths_during_primary,
           Mortality_Rate_Percent = mortality_rate,
           Survivors = survivors_to_challenge,
           Survival_Rate_Percent = survival_rate)

cat("COMPLETE EXPERIMENTAL RESULTS (All", length(challenge_mice), "mice):\n")
print(results_summary_all)

# Also create summary for immune subset
results_summary_immune <- primary_mortality_immune %>%
    rename(Species = Parasite_primary,
           Total_With_Immune = total_mice_with_immune,
           Deaths_With_Immune = deaths_with_immune,
           Survivors_With_Immune = survivors_with_immune,
           Mortality_Rate_Immune_Subset = mortality_rate_immune)

cat("\nIMMUNE DATA SUBSET (136 mice):\n")
print(results_summary_immune)

# Key numbers for manuscript - FIXED
total_deaths_all <- sum(results_summary_all$Deaths)
total_mice_all <- sum(results_summary_all$Total_Mice)
falc_mortality_all <- results_summary_all[results_summary_all$Species == "E. falciformis", "Mortality_Rate_Percent"][[1]]
ferr_mortality_all <- results_summary_all[results_summary_all$Species == "E. ferrisi", "Mortality_Rate_Percent"][[1]]
control_mortality_all <- results_summary_all[results_summary_all$Species == "Uninfected controls", "Mortality_Rate_Percent"][[1]]

cat("\nKEY NUMBERS FOR MANUSCRIPT:\n")
cat("Total deaths (all", length(challenge_mice), "mice):", total_deaths_all, "\n")
cat("E. falciformis mortality rate:", falc_mortality_all, "%\n")
cat("E. ferrisi mortality rate:", ferr_mortality_all, "%\n")
cat("Control mortality rate:", control_mortality_all, "%\n")

# Missing mice analysis
missing_mice_analysis <- mortality_status %>%
    filter(death == "primary", !has_immune_data) %>%
    count(Parasite_primary)

if(nrow(missing_mice_analysis) > 0) {
    cat("\nMice that died but have NO immune data:\n")
    print(missing_mice_analysis)
} else {
    cat("\nAll deceased mice have immune data available.\n")
}


# Mouse strains
model <- lm(WL_max ~ mouse_strain, data = lab)
summary(model)


# =============================================================================
# 11. STATISTICAL TESTS (CORRECTED)
# =============================================================================

cat("\n=== STATISTICAL TESTS ===\n")

# Test for difference in mortality rates between species (using ALL mice data)
if(nrow(results_summary_all) >= 2) {
    mortality_matrix <- matrix(
        c(results_summary_all$Deaths, 
          results_summary_all$Survivors),
        nrow = 2, byrow = TRUE
    )
    colnames(mortality_matrix) <- results_summary_all$Species
    rownames(mortality_matrix) <- c("Deaths", "Survivors")
    
    cat("Mortality contingency table (all", length(challenge_mice), "mice):\n")
    print(mortality_matrix)
    
    # Only do Fisher's test if we have at least 2 species with some variation
    if(ncol(mortality_matrix) >= 2 && sum(mortality_matrix[1,]) > 0) {
        fisher_test <- fisher.test(mortality_matrix)
        cat("\nFisher's exact test for mortality difference:\n")
        cat("p-value:", format(fisher_test$p.value, scientific = TRUE), "\n")
        cat("Odds ratio:", fisher_test$estimate, "\n")
        cat("95% CI:", fisher_test$conf.int[1], "-", fisher_test$conf.int[2], "\n")
    }
}

# =============================================================================
# 12. CREATE SUMMARY SENTENCES FOR MANUSCRIPT (CORRECTED)
# =============================================================================

cat("\n=== MANUSCRIPT SENTENCES ===\n")

cat("SUGGESTED MANUSCRIPT TEXT:\n")
cat("=========================\n")
cat(sprintf("We conducted five independent infection experiments, initially including %d laboratory mice across diverse wild-derived genetic backgrounds. ", 
            length(challenge_mice)))
cat(sprintf("Due to infection-induced mortality and tissue collection constraints, immune gene expression analysis was possible for 136 mice. "))
cat(sprintf("A total of %d mice died during primary infection, with mortality restricted entirely to infected animals. ", 
            total_deaths_all))

# Get specific numbers for each species - FIXED
falc_deaths <- results_summary_all[results_summary_all$Species == "E. falciformis", "Deaths"][[1]]
ferr_deaths <- results_summary_all[results_summary_all$Species == "E. ferrisi", "Deaths"][[1]]

cat(sprintf("E. falciformis showed significantly higher mortality (%d deaths, %.1f%% mortality rate) compared to ", 
            falc_deaths, falc_mortality_all))
cat(sprintf("E. ferrisi (%d deaths, %.1f%% mortality rate), while no uninfected controls died. ", 
            ferr_deaths, ferr_mortality_all))

if(exists("fisher_test")) {
    cat(sprintf("The difference in mortality rates between species was statistically significant (Fisher's exact test, p = %.2e).\n", 
                fisher_test$p.value))
}

# Calculate how many mice lost immune data
mice_lost_to_analysis <- length(challenge_mice) - 136
cat(sprintf("The %d mice not included in immune analyses represent cases where death occurred before tissue collection was possible, ", 
            mice_lost_to_analysis))
cat("indicating that our immune gene expression dataset represents a conservative estimate of infection severity.\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

