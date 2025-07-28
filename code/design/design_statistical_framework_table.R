# =============================================
# TABLE: STATISTICAL ANALYSIS FRAMEWORK
# =============================================

# Create the statistical analysis framework table
create_statistical_framework_table <- function() {
    
    # Define the table data
    statistical_framework <- tribble(
        ~Phase, ~Research_Question, ~Analysis, ~Method, ~Key_Result, ~Sample_Size, ~Performance,
        
        # Laboratory Development Phase
        "Discovery", "Can immune genes predict infection costs?", "DISC-1", "Linear regression (PC1, PC2 → weight loss)", "Significant but modest prediction", "n = 136", "R² = 0.106***",
        "Optimization", "Can machine learning improve prediction?", "DISC-2", "Random forest (19 genes → weight loss)", "Substantial improvement achieved", "n = 136", "R² = 0.476***",
        "Validation", "Is the model reliable?", "DISC-3", "Train-test cross-validation", "Strong predictive accuracy", "n = 95→41", "r = 0.79***",
        
        # Field Translation Phase  
        "Detection", "Does the model work in wild populations?", "FIELD-1", "Predicted vs. observed infection status", "Successfully detects infection", "n = 305", "+1.15%***",
        "Discrimination", "Can it distinguish parasite species?", "FIELD-2", "Predicted loss by species identity", "Species-specific responses", "n = 169", "E.f: +2.06%**, E.r: +1.25%**",
        "Scaling", "Does it correlate with infection severity?", "FIELD-3", "Predicted loss vs. parasite load", "Scales with infection intensity", "n = 185", "r = 0.233*",
        
        # Biological Validation Phase
        "Physiological relevance", "Does it capture real health impacts?", "PROOF-1", "Predicted loss vs. body condition", "Correlates with actual body weight", "n = 336", "ρ = -0.115*",
        "Specificity", "Is the response Eimeria-specific?", "PROOF-2", "Predicted loss vs. parasite community", "Specific to Eimeria infections only", "n = 305", "p < 0.001***"
    ) %>%
        mutate(
            # Add phase groupings
            Phase_Group = case_when(
                Phase %in% c("Discovery", "Optimization", "Validation") ~ "Laboratory Development",
                Phase %in% c("Detection", "Discrimination", "Scaling") ~ "Field Translation", 
                Phase %in% c("Physiological relevance", "Specificity") ~ "Biological Validation"
            )
        ) %>%
        # Reorder columns
        dplyr::select(Phase_Group, Phase, Research_Question, Analysis, Method, Key_Result, Sample_Size, Performance)
    
    return(statistical_framework)
}

# Create the table
statistical_framework_data <- create_statistical_framework_table()

# Create the publication-ready gt table
statistical_framework_table <- statistical_framework_data %>%
    gt(groupname_col = "Phase_Group") %>%
    tab_header(
        title = "",
        subtitle = ""
    ) %>%
    cols_label(
        Phase = "Analysis Phase",
        Research_Question = "Research Question", 
        Analysis = "Model ID",
        Method = "Statistical Method",
        Key_Result = "Key Finding",
        Sample_Size = "Sample Size",
        Performance = "Performance Metric"
    ) %>%
    # Style the group headers (phase groups)
    tab_style(
        style = list(
            cell_fill(color = "#E8F4FD"),
            cell_text(weight = "bold", size = px(14))
        ),
        locations = cells_row_groups()
    ) %>%
    # Style column headers
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Style the Model ID column
    tab_style(
        style = list(
            cell_text(weight = "bold", color = "#1f77b4"),
            cell_fill(color = "#f8f9fa")
        ),
        locations = cells_body(columns = Analysis)
    ) %>%
    # Style performance metrics with significance
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = Performance)
    ) %>%
    # Add footnotes
    tab_footnote(
        footnote = "Significance levels: *p < 0.05, **p < 0.01, ***p < 0.001",
        locations = cells_column_labels(columns = Performance)
    ) %>%
    tab_footnote(
        footnote = "E.f: Eimeria falciformis; E.r: E. ferrisi",
        locations = cells_body(columns = Performance, rows = 5)
    ) %>%
    tab_footnote(
        footnote = "Progressive sample sizes reflect train→test validation approach",
        locations = cells_body(columns = Sample_Size, rows = 3)
    ) %>%
    # Adjust column widths
    cols_width(
        Phase ~ px(120),
        Research_Question ~ px(200),
        Analysis ~ px(80),
        Method ~ px(250),
        Key_Result ~ px(200),
        Sample_Size ~ px(80),
        Performance ~ px(120)
    ) %>%
    # Style significance asterisks
    tab_style(
        style = cell_text(color = "#d62728", weight = "bold"),
        locations = cells_body(
            columns = Performance,
            rows = str_detect(Performance, "\\*")
        )
    ) %>%
    # Add source note
    tab_source_note(
        source_note = "Framework progresses from basic linear prediction (R² = 0.106) through machine learning optimization (R² = 0.476) to comprehensive field validation with biological relevance."
    )

# Print the table
print(statistical_framework_table)

# Save the table using your existing function
save_table_all_formats(statistical_framework_table, "Table_Statistical_Analysis_Framework")

# Alternative version with simpler formatting for supplementary material
statistical_framework_simple <- statistical_framework_data %>%
    gt() %>%
    tab_header(
        title = "Table 1. Statistical Analysis Framework",
        subtitle = "Complete analytical approach from laboratory development to field validation"
    ) %>%
    cols_label(
        Phase_Group = "Phase",
        Phase = "Step",
        Research_Question = "Research Question",
        Analysis = "Model",
        Method = "Method", 
        Key_Result = "Result",
        Sample_Size = "n",
        Performance = "Performance"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    fmt_markdown(columns = c(Method, Key_Result)) %>%
    tab_footnote(
        footnote = "*p < 0.05, **p < 0.01, ***p < 0.001",
        locations = cells_column_labels(columns = Performance)
    )

print(statistical_framework_simple)

# Save simple version too
save_table_all_formats(statistical_framework_simple, "Table_Statistical_Framework_Simple")

# Create a summary version for methods section
methods_summary <- statistical_framework_data %>%
    group_by(Phase_Group) %>%
    summarise(
        Models = paste(unique(Analysis), collapse = ", "),
        Primary_Questions = n(),
        Total_Sample_Sizes = paste(unique(Sample_Size), collapse = "; "),
        .groups = "drop"
    ) %>%
    gt() %>%
    tab_header(
        title = "Statistical Analysis Summary by Phase",
        subtitle = "Overview of analytical framework for methods section reference"
    ) %>%
    cols_label(
        Phase_Group = "Analysis Phase",
        Models = "Model IDs", 
        Primary_Questions = "Number of Tests",
        Total_Sample_Sizes = "Sample Sizes"
    )

print(methods_summary)

save_table_all_formats(methods_summary, "Table_Methods_Summary")

# Print confirmation
cat("✅ Statistical Analysis Framework tables created and saved!\n")
cat("   - Main table: Table_Statistical_Analysis_Framework\n") 
cat("   - Simple version: Table_Statistical_Framework_Simple\n")
cat("   - Methods summary: Table_Methods_Summary\n")
cat("   - All saved in multiple formats (HTML, DOCX, PNG, PDF, TEX)\n")

