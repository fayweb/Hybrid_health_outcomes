# ***********************************************************
# design_field_quantification_methods.R
# Purpose: Summarize and document how Eimeria infection status
# was determined in field mice for downstream analyses.
# ***********************************************************
# ---- Import cleaned data ----
hm <- read.csv("data/analysis/final/hm_ready_for_analysis.csv")


# Filter to field samples only
field <- hm %>%
    filter(origin == "Field")

cat("=== FIELD DATA OVERVIEW ===\n")
cat("Total field mice with immune data:", length(unique(field$Mouse_ID)), "\n\n")

# Define gene vector from main pipeline
Genes_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
             "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
             "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
             "TICAM1", "TNF")

# Check availability
available_genes <- intersect(Genes_v, colnames(field))
cat("Available immune genes:", length(available_genes), "out of", length(Genes_v), "\n")
cat("Missing genes:", setdiff(Genes_v, available_genes), "\n\n")

# Deduplicate and flag immune data availability
field <- field %>%
    group_by(Mouse_ID) %>%
    summarise(across(everything(), first), .groups = "drop") %>%
    mutate(has_immune_data = TRUE)

cat("Unique mice after deduplication:", nrow(field), "\n")

# ---- Check data availability for different methods ----
field <- field %>%
    mutate(
        has_qPCR_tissue = !is.na(delta_ct_cewe_MminusE),
        has_qPCR_feces = !is.na(FEC_Eim_Ct),
        has_oocyst_counts = !is.na(OPG) & OPG >= 0,
        has_amplicon_data = !is.na(amplicon_species),
        has_species_ID = !is.na(species_Eimeria),  # Fixed: removed != "NA"
        has_infection_status = !is.na(MC.Eimeria)  # Added for clarity
    )

# Summary of data availability
data_availability <- field %>%
    summarise(
        total_mice = n(),
        has_immune_genes = sum(has_immune_data),
        has_qPCR_tissue = sum(has_qPCR_tissue),
        has_qPCR_feces = sum(has_qPCR_feces),
        has_oocyst_counts = sum(has_oocyst_counts),
        has_amplicon_data = sum(has_amplicon_data),
        has_species_ID = sum(has_species_ID),
        has_infection_status = sum(has_infection_status)
    )

print(data_availability)

# ---- Derive infection status and species following hierarchical logic ----
field <- field %>%
    mutate(
        # Infection status from qPCR tissue
        infection_tissue = case_when(
            MC.Eimeria == "TRUE" ~ "TRUE",
            MC.Eimeria == "FALSE" ~ "FALSE",
            TRUE ~ "Unknown"
        ),
        
        # Infection status from species identification
        infection_species = case_when(
            species_Eimeria %in% c("E. ferrisi", "E. falciformis") ~ "TRUE",
            species_Eimeria == "Uninfected" ~ "FALSE",
            TRUE ~ "Unknown"
        ),
        
        # Final infection status (prioritize qPCR, then species)
        infection_status = case_when(
            infection_tissue %in% c("TRUE", "FALSE") ~ infection_tissue,
            infection_species %in% c("TRUE", "FALSE") ~ infection_species,
            TRUE ~ "Unknown"
        ),
        
        # Clean species names
        species_clean = case_when(
            species_Eimeria == "E. ferrisi" ~ "E. ferrisi",
            species_Eimeria == "E. falciformis" ~ "E. falciformis", 
            species_Eimeria == "Uninfected" ~ "Uninfected",
            TRUE ~ "Unknown"
        )
    )

# ---- Identify method used for infection status assignment ----
field <- field %>%
    mutate(status_source = case_when(
        !is.na(MC.Eimeria) ~ "qPCR (tissue)",
        is.na(MC.Eimeria) & !is.na(species_Eimeria) ~ "Species call (qPCR/amplicon)",
        TRUE ~ "Unknown"
    ))

# ---- Define mice usable for downstream modeling ----
field <- field %>%
    mutate(usable_for_model = infection_status %in% c("TRUE", "FALSE") &
               species_clean %in% c("E. ferrisi", "E. falciformis", "Uninfected"))

cat("\n=== Mice usable for downstream models ===\n")
table(field$usable_for_model)

# ---- Create supplementary tables ----

# Supplementary Table S1: Methods summary
methods_used <- tibble::tibble(
    Method = c("Tissue qPCR (detection)", "Tissue qPCR (species)", "Amplicon sequencing"),
    Variable = c("MC.Eimeria", "eimeriaSpecies", "amplicon_species"),
    Purpose = c(
        "Primary infection detection (presence/absence)",
        "Primary species identification",
        "Backup species identification when qPCR fails"
    ),
    Final_Variable = c("infection_status", "species_Eimeria", "species_Eimeria"),
    Mice_with_Data = c(
        sum(!is.na(field$MC.Eimeria)),
        sum(!is.na(field$species_Eimeria)),
        sum(!is.na(field$amplicon_species))
    )
)

cat("\n=== Supplementary Table S1: Methods Summary ===\n")
print(methods_used)

# Supplementary Table S2: Final infection status
infection_summary <- field %>%
    filter(usable_for_model == TRUE) %>%
    count(infection_status) %>%
    rename(Number_of_Mice = n)

cat("\n=== Supplementary Table S2: Final Infection Status ===\n")
print(infection_summary)

# Supplementary Table S3: Species identification
species_summary <- field %>%
    filter(usable_for_model == TRUE) %>%
    count(species_clean) %>%
    rename(Number_of_Mice = n)

cat("\n=== Supplementary Table S3: Species Identification ===\n")
print(species_summary)

# Supplementary Table S4: Method source breakdown
method_source_summary <- field %>%
    filter(usable_for_model == TRUE) %>%
    count(status_source) %>%
    rename(Number_of_Mice = n)

cat("\n=== Supplementary Table S4: Data Source for Final Assignment ===\n")
print(method_source_summary)

# ---- Diagnostic checks ----
cat("\n=== DIAGNOSTIC CHECKS ===\n")
cat("MC.Eimeria distribution:\n")
print(table(field$MC.Eimeria, useNA = "always"))

cat("\nspecies_Eimeria distribution:\n")
print(table(field$species_Eimeria, useNA = "always"))

cat("\nFinal infection_status distribution:\n")
print(table(field$infection_status, useNA = "always"))

cat("\nStatus source breakdown:\n")
print(table(field$status_source, useNA = "always"))


# =============================================
# SUPPLEMENTARY TABLE: EIMERIA DETECTION METHODS
# =============================================
# Create methods summary table
methods_table_data <- tibble(
    Method = c(
        "Caecal qPCR + melting curve", 
        "Caecal qPCR + melting curve", 
        "Amplicon sequencing"
    ),
    Variable = c("MC.Eimeria", "eimeriaSpecies", "amplicon_species"),
    Purpose = c(
        "Infection detection (presence/absence)",
        "Species identification", 
        "Backup species identification"
    ),
    Priority = c("Primary", "Primary", "Backup"),
    Mice_with_Data = c(185, 169, 134),
    Total_Mice = c(336, 336, 336),
    Success_Rate = paste0(round((c(185, 169, 134) / 336) * 100, 1), "%")
)

supp_table_methods <- methods_table_data %>%
    gt() %>%
    tab_header(
        title = "Eimeria detection methods in field mice",
        subtitle = "Summary of hierarchical approach for infection status and species assignment"
    ) %>%
    cols_label(
        Method = "Detection Method",
        Variable = "Variable Name",
        Purpose = "Purpose",
        Priority = "Priority",
        Mice_with_Data = "Mice with Data",
        Total_Mice = "Total Mice",
        Success_Rate = "Success Rate"
    ) %>%
    # Combine sample size columns
    cols_merge(
        columns = c(Mice_with_Data, Total_Mice),
        pattern = "{1}/{2}"
    ) %>%
    cols_label(
        Mice_with_Data = "Sample Size (n/total)"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Highlight primary vs backup methods
    tab_style(
        style = cell_fill(color = "#f3e5f5"),
        locations = cells_body(rows = Priority == "Primary")
    ) %>%
    tab_style(
        style = cell_fill(color = "#fff3e0"), 
        locations = cells_body(rows = Priority == "Backup")
    ) %>%
    tab_footnote(
        footnote = "Primary method: Caecal tissue qPCR with melting curve analysis for both detection and species ID",
        locations = cells_column_labels(columns = Method)
    ) %>%
    tab_footnote(
        footnote = "Backup method: Used only when qPCR species identification failed",
        locations = cells_body(columns = Method, rows = 3)
    )

print(supp_table_methods)

# =============================================
# SUPPLEMENTARY TABLE: FINAL SAMPLE COMPOSITION
# =============================================

# Sample composition data
composition_data <- tibble(
    Category = c(
        "Total field mice", 
        "Usable for modeling",
        "Excluded from modeling",
        "",  # spacer
        "Infected",
        "Uninfected", 
        "",  # spacer
        "E. ferrisi",
        "E. falciformis", 
        "Uninfected"
    ),
    Number_of_Mice = c(336, 169, 167, NA, 65, 104, NA, 45, 14, 110),
    Percentage = c(
        "100%", "50.3%", "49.7%", "", 
        "38.5%*", "61.5%*", "",
        "26.6%*", "8.3%*", "65.1%*"
    ),
    Group = c(
        "Overview", "Overview", "Overview", "",
        "Infection Status", "Infection Status", "",
        "Species Identity", "Species Identity", "Species Identity"
    )
)

supp_table_composition <- composition_data %>%
    filter(!is.na(Number_of_Mice)) %>%  # Remove spacer rows
    gt(groupname_col = "Group") %>%
    tab_header(
        title = "Supplementary Table S2. Final sample composition for analyses",
        subtitle = "Breakdown of 336 field mice by infection status and species"
    ) %>%
    cols_label(
        Category = "Category",
        Number_of_Mice = "Number of Mice", 
        Percentage = "Percentage"
    ) %>%
    # Format species names in italics
    text_transform(
        locations = cells_body(columns = Category),
        fn = function(x) {
            case_when(
                str_detect(x, "E\\. ") ~ paste0("<em>", x, "</em>"),
                TRUE ~ x
            )
        }
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
    ) %>%
    tab_footnote(
        footnote = "* Percentages calculated from usable mice (n=169)",
        locations = cells_body(columns = Percentage, rows = str_detect(Percentage, "\\*"))
    ) %>%
    fmt_markdown(columns = Category)

print(supp_table_composition)

# =============================================
# SUPPLEMENTARY TABLE: DATA SOURCE BREAKDOWN
# =============================================

source_data <- tibble(
    Assignment_Method = c("qPCR (tissue)", "Species call (qPCR/amplicon)"),
    Number_of_Mice = c(49, 120),
    Percentage = c("29.0%", "71.0%"),
    Description = c(
        "Direct infection status from MC.Eimeria",
        "Infection status inferred from species identification"
    )
)

supp_table_sources <- source_data %>%
    gt() %>%
    tab_header(
        title = "Supplementary Table S3. Data source for final infection assignments",
        subtitle = "Method used to determine infection status in usable mice (n=169)"
    ) %>%
    cols_label(
        Assignment_Method = "Assignment Method",
        Number_of_Mice = "Number of Mice",
        Percentage = "Percentage",
        Description = "Description"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    tab_footnote(
        footnote = "Species identification prioritized qPCR melting curves, with amplicon sequencing as backup",
        locations = cells_body(columns = Assignment_Method, rows = 2)
    )

print(supp_table_sources)

# =============================================
# SAVE ALL TABLES
# =============================================

# Save using your existing function
save_table_all_formats(supp_table_methods, "Supplementary_Table_S1_Eimeria_methods")
save_table_all_formats(supp_table_composition, "Supplementary_Table_S2_sample_composition") 
save_table_all_formats(supp_table_sources, "Supplementary_Table_S3_data_sources")

cat("✅ All supplementary tables saved successfully!\n")
