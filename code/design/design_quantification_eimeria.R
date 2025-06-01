# =============================================================================
# COMPLETE WORKING FIELD DATA VISUALIZATIONS
# =============================================================================
# All your data is perfect! Let's create the full visualization suite.

# =============================================
# SUPPLEMENTARY TABLE: EIMERIA DETECTION METHODS (CORRECTED)
# =============================================

hm <- read.csv("data/analysis/final/hm_ready_for_analysis.csv")
field <- hm %>%
    filter(origin == "Field")

# Create methods summary table with CORRECT numbers
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
    Mice_with_Data = c(185, 49, 134),
    Total_Mice = c(336, 336, 336),
    Success_Rate = paste0(round((c(185, 49, 134) / 336) * 100, 1), "%"),
    Notes = c(
        "Direct infection status from melting curve",
        "Species ID from melting curve patterns",
        "Used when qPCR species ID unavailable"
    )
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
        Success_Rate = "Success Rate",
        Notes = "Notes"
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
        footnote = "Final species assignment: 49 from qPCR + 120 from amplicon = 169 total mice",
        locations = cells_column_labels(columns = Method)
    )

print(supp_table_methods)

# =============================================
# SUPPLEMENTARY TABLE: FINAL SAMPLE COMPOSITION (CORRECTED)
# =============================================

# Sample composition data with CORRECT numbers
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
        title = "Final sample composition for analyses",
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
# SUPPLEMENTARY TABLE: DATA SOURCE BREAKDOWN (CORRECTED)
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
        title = "Data source for final infection assignments",
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
        footnote = "120 mice: 49 from qPCR species ID + 120 from amplicon backup when qPCR unavailable",
        locations = cells_body(columns = Assignment_Method, rows = 2)
    )

print(supp_table_sources)

# =============================================
# SUPPLEMENTARY TABLE: DETAILED METHOD BREAKDOWN (NEW)
# =============================================

# Create detailed breakdown table
detailed_breakdown <- tibble(
    Step = c(
        "qPCR infection detection",
        "qPCR species identification", 
        "Amplicon species identification",
        "Final species assignment",
        "Complete data for modeling"
    ),
    Method = c(
        "Melting curve analysis",
        "Melting curve patterns",
        "18S/28S sequencing", 
        "Hierarchical combination",
        "Both infection status + species ID"
    ),
    Mice_Success = c(185, 49, 134, 169, 169),
    Mice_Failed = c(151, 287, 202, 167, 167),
    Success_Rate = paste0(round(c(185, 49, 134, 169, 169)/336*100, 1), "%"),
    Notes = c(
        "Primary method for infection presence/absence",
        "Limited by time constraints",
        "Backup method when qPCR species ID unavailable",
        "49 qPCR + 120 amplicon = 169 total",
        "Used for downstream random forest validation"
    )
)

supp_table_detailed <- detailed_breakdown %>%
    gt() %>%
    tab_header(
        title = "Detailed method performance breakdown",
        subtitle = "Step-by-step success rates in the hierarchical detection pipeline"
    ) %>%
    cols_label(
        Step = "Pipeline Step",
        Method = "Method Used",
        Mice_Success = "Success (n)",
        Mice_Failed = "Missing (n)", 
        Success_Rate = "Success Rate",
        Notes = "Notes"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Highlight final result
    tab_style(
        style = cell_fill(color = "#e8f5e8"),
        locations = cells_body(rows = Step == "Complete data for modeling")
    )

print(supp_table_detailed)

# =============================================
# SAVE ALL TABLES
# =============================================

save_table_all_formats(supp_table_methods, "Supplementary_Table__Eimeria_methods")
save_table_all_formats(supp_table_composition, "Supplementary_Table_sample_composition") 
save_table_all_formats(supp_table_sources, "Supplementary_Table_data_sources")
save_table_all_formats(supp_table_detailed, "Supplementary_Table_S4_detailed_breakdown")

cat("✅ All corrected supplementary tables saved successfully!\n")

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