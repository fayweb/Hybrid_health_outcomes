
# =============================================
# IMPROVED TABLE: STATISTICAL ANALYSIS FRAMEWORK
# Including TRANS-1 and corrected workflow
# =============================================

library(gt)
library(tibble)
library(dplyr)
library(stringr)

# Create the improved statistical analysis framework table
create_improved_statistical_framework_table <- function() {
    
    # Define the table data with complete workflow
    statistical_framework <- tribble(
        ~Phase, ~Research_Question, ~Analysis, ~Method, ~Key_Result, ~Sample_Size, ~Performance,
        
        # Laboratory Development Phase
        "Discovery", "Can immune genes predict infection costs?", "DISC-1", "Linear regression (PC1, PC2 → weight loss)", "Significant but modest prediction", "n = 136", "R² = 0.106***",
        "Optimization", "Can machine learning improve prediction?", "DISC-2", "Random forest (19 genes → weight loss)", "Substantial improvement achieved", "n = 136", "R² = 0.476***",
        "Validation", "Is the model reliable?", "DISC-3", "Train-test cross-validation", "Strong predictive accuracy", "n = 136 (70/30 split)", "r = 0.79***",
        
        # Cross-Population Translation Phase (NEW!)
        "Gene Validation", "Which genes show consistent responses across populations?", "TRANS-1", "Linear regression per gene (lab vs field)", "3 genes cross-validated", "n = 305", "3/19 genes validated",
        
        # Field Translation Phase  
        "Detection", "Does the model work in wild populations?", "FIELD-1", "Predicted vs. observed infection status", "Successfully detects infection", "n = 305", "+1.15%***",
        "Discrimination", "Can it distinguish parasite species?", "FIELD-2", "Predicted loss by species identity", "Species-specific responses", "n = 169", "E.f: +2.06%**, E.r: +1.25%**",
        "Scaling", "Does it correlate with infection severity?", "FIELD-3", "Predicted loss vs. parasite load", "Scales with infection intensity", "n = 185", "r = 0.233*",
        
        # Biological Validation Phase
        "Physiological relevance", "Does it capture real health impacts?", "PROOF-1", "Predicted loss vs. body condition", "Correlates with actual body weight", "n = 336", "ρ = -0.115*",
        "Specificity", "Is the response Eimeria-specific?", "PROOF-2", "Predicted loss vs. parasite community", "Specific to Eimeria infections only", "n = 305", "p < 0.001***"
    ) %>%
        mutate(
            # Add phase groupings with improved logic
            Phase_Group = case_when(
                Phase %in% c("Discovery", "Optimization", "Validation") ~ "Laboratory Development",
                Phase %in% c("Gene Validation") ~ "Cross-Population Translation",
                Phase %in% c("Detection", "Discrimination", "Scaling") ~ "Field Translation", 
                Phase %in% c("Physiological relevance", "Specificity") ~ "Biological Validation"
            )
        ) %>%
        # Reorder columns
        dplyr::select(Phase_Group, Phase, Research_Question, Analysis, Method, Key_Result, Sample_Size, Performance)
    
    return(statistical_framework)
}

# Create the table
statistical_framework_data <- create_improved_statistical_framework_table()

# Create the publication-ready gt table with improved styling
statistical_framework_table <- statistical_framework_data %>%
    gt(groupname_col = "Phase_Group") %>%
    tab_header(
        title = "Table 1. Statistical Analysis Framework",
        subtitle = "Complete analytical workflow from laboratory development through field validation to biological proof-of-concept"
    ) %>%
    cols_label(
        Phase = "Analysis Phase",
        Research_Question = "Research Question", 
        Analysis = "Model ID",
        Method = "Statistical Method",
        Key_Result = "Key Finding",
        Sample_Size = "Sample Size",
        Performance = "Performance Metric¹"
    ) %>%
    # Enhanced styling for group headers
    tab_style(
        style = list(
            cell_fill(color = "#E3F2FD"),
            cell_text(weight = "bold", size = px(13), color = "#1565C0")
        ),
        locations = cells_row_groups()
    ) %>%
    # Style column headers
    tab_style(
        style = list(
            cell_text(weight = "bold", size = px(12)),
            cell_fill(color = "#F5F5F5")
        ),
        locations = cells_column_labels()
    ) %>%
    # Highlight Model ID column
    tab_style(
        style = list(
            cell_text(weight = "bold", color = "#1976D2"),
            cell_fill(color = "#FAFAFA")
        ),
        locations = cells_body(columns = Analysis)
    ) %>%
    # Style performance metrics with significance
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = Performance)
    ) %>%
    # Highlight significant results
    tab_style(
        style = cell_text(color = "#D32F2F", weight = "bold"),
        locations = cells_body(
            columns = Performance,
            rows = str_detect(Performance, "\\*")
        )
    ) %>%
    # Add comprehensive footnotes
    tab_footnote(
        footnote = "Significance levels: *p < 0.05, **p < 0.01, ***p < 0.001",
        locations = cells_column_labels(columns = Performance)
    ) %>%
    tab_footnote(
        footnote = "E.f: Eimeria falciformis; E.r: E. ferrisi",
        locations = cells_body(columns = Performance, rows = 6)
    ) %>%
    tab_footnote(
        footnote = "Train-test validation used 70% training, 30% testing from full dataset",
        locations = cells_body(columns = Sample_Size, rows = 3)
    ) %>%
    tab_footnote(
        footnote = "Cross-validated genes: CXCL9 (both species), TICAM1, PRF1 (E. falciformis)",
        locations = cells_body(columns = Key_Result, rows = 4)
    ) %>%
    # Adjust column widths for better readability
    cols_width(
        Phase ~ px(140),
        Research_Question ~ px(220),
        Analysis ~ px(90),
        Method ~ px(280),
        Key_Result ~ px(220),
        Sample_Size ~ px(100),
        Performance ~ px(140)
    ) %>%
    # Add comprehensive source note
    tab_source_note(
        source_note = "Framework demonstrates progression from basic linear prediction (R² = 0.106) through machine learning optimization (R² = 0.476) to comprehensive field validation with biological relevance. Cross-population translation validates 3/19 genes as conserved biomarkers."
    ) %>%
    # Format text properly
    fmt_markdown(columns = c(Method, Key_Result))

# Print the table
print(statistical_framework_table)

# Create a methods-focused summary table
methods_focused_summary <- statistical_framework_data %>%
    group_by(Phase_Group) %>%
    summarise(
        Models_Included = paste(unique(Analysis), collapse = ", "),
        Key_Questions = n(),
        Sample_Range = paste(range(parse_number(Sample_Size)), collapse = "-"),
        Primary_Outcome = case_when(
            Phase_Group == "Laboratory Development" ~ "Model Development & Validation",
            Phase_Group == "Cross-Population Translation" ~ "Gene Conservation Assessment", 
            Phase_Group == "Field Translation" ~ "Field Application Success",
            Phase_Group == "Biological Validation" ~ "Biological Relevance Proof"
        ),
        .groups = "drop"
    ) %>%
    gt() %>%
    tab_header(
        title = "Statistical Analysis Summary by Phase",
        subtitle = "Methodological overview for manuscript methods section"
    ) %>%
    cols_label(
        Phase_Group = "Analysis Phase",
        Models_Included = "Model IDs", 
        Key_Questions = "Number of Tests",
        Sample_Range = "Sample Size Range",
        Primary_Outcome = "Phase Objective"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    )

# Print methods summary
print(methods_focused_summary)

# Save both tables using your existing function (assuming it exists)
 save_table_all_formats(statistical_framework_table, "Table_1_Statistical_Framework_Complete")
 save_table_all_formats(methods_focused_summary, "Table_Methods_Summary_Improved")


