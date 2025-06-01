# ***********************************************************
# Title: Predicting the health outcomes of infections in hybrid mice

# Purpose: This script analysis the potential of each immune 
# gene used int the experiment to predict the 
# the weight loss during infection
# Authors: Fay Webster


# WOrking with laboratory data only
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

# List of dependent variables
dependent_vars <- c("IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10",
                    "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                    "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                    "TICAM1", "TNF")

# Assuming 'lab' is your dataframe
#lab$current_infection <- factor(lab$current_infection, levels = c("uninfected", "E_falciformis", "E_ferrisi"))

# Perform regressions
results <- lapply(dependent_vars, function(var) {
    lm_formula <- as.formula(paste(var, "~ current_infection", sep = " "))
    lm(lm_formula, data = lab)
})

names(results) <- dependent_vars  # Name each regression with the name of the dependent variable


# Correctly add 'model' column to each tidied model's dataframe
tidy_models <- do.call(rbind, lapply(names(results), function(name) {
    model_df <- tidy(results[[name]], conf.int = TRUE)
    model_df$model <- name  # Add model name as a new column
    return(model_df)
}))

# Now 'tidy_models' should have a 'model' column correctly populated

# Filter out intercept terms
tidy_models_no_intercept <- tidy_models[!grepl("intercept", tidy_models$term, ignore.case = TRUE),]

tidy_models_no_intercept <- as.data.frame(tidy_models_no_intercept) %>%
    mutate(term = case_when(
        term == "current_infectionE. falciformis" ~ "E. falciformis",
        term == "current_infectionE. ferrisi" ~ "E. ferrisi",
    ))




# Determine the common y-axis range
common_y_limits <- range(-7,7)

coef_mmr <- ggplot(tidy_models_no_intercept, aes(x = model, y = estimate, color = term)) +
    coord_flip() +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("E. falciformis" = "#FF0000", "E. ferrisi" = "#7A0092")) +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Coefficients estimate (Difference to uninfected)") +
    scale_y_continuous(limits = common_y_limits) +  # Apply shared y-axis limits
    theme(legend.title = element_blank(),
          legend.position = "none")


print(coef_mmr)

ggsave(filename = paste0(an_fi, "/coef_plot_lab_genes.pdf"),
       plot = coef_mmr, width = 6, height = 4, dpi = 300)



####################################################
####################################################
############################# gene expression distribution
lab %>%
    pivot_longer(cols = all_of(dependent_vars), 
                 names_to = "Genes", values_to = "Expression") %>%
    ggplot(aes(x = Expression, fill = current_infection)) +
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) +
    geom_boxplot(position = "dodge2",
                 width = .5, 
                 outlier.shape = NA,
                 orientation = "y"  # Set orientation to y
    ) +
    facet_wrap(~Genes,  scales = 'free', ncol = 4) +
    labs(x = "Expression Level", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = color_mapping, labels = labels)  +
    theme(legend.title = element_blank(), 
          legend.position = c(0.85, 0.06),
          legend.text = element_markdown())+
    labs(y = "Density", 
         x = "Gene expression level") -> density_imm

density_imm

ggsave(filename = paste0(an_fi, "/density_immune_genes.jpeg"),
       plot = density_imm, width = 10, height = 8, dpi = 300)


#######################
# combine
comb <- (density_imm | coef_mmr) +
    #  plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
comb <- comb + 
    plot_annotation(title = 'Fig. 2', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))

comb <- comb + 
    plot_layout(heights = c(1, 1), 
                widths = c(2, 1))



# Display the panel figure
# print(comb)

# Save the panel figure
ggsave(paste0(panels_fi, "/panel_immune_gene_expression_lab.jpeg"), 
       comb, width = 12, height = 6, dpi = 300)    



# Update gene names to replace "." with "-" and "IFNg" with "IFNy"
regression_table <- tidy_models_no_intercept %>%
    dplyr::select(model, term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(
        model = gsub("\\.", "-", model),  # Replace "." with "-"
        model = gsub("IFNg", "IFNy", model),  # Replace "IFNg" with "IFNy"
        significance = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
        )
    ) %>%
    rename(
        Gene = model,
        'Treatment group' = term,
        Estimate = estimate,
        `Std. Error` = std.error,
        `CI Lower` = conf.low,
        `CI Upper` = conf.high,
        `P-value` = p.value,
        `Significance` = significance
    )

caption_text <- "Table 1: Regression results for immune gene expression in E. ferrisi and E. falciformis-infected mice in controlled experiments.<br><em>Significance codes:</em> '***' 0.001, '**' 0.01, '*' 0.05"

publication_table <- regression_table %>%
    kbl(
        caption = caption_text,
        format = "html",
        escape = FALSE,
        digits = 3
    ) %>%
    kable_styling(
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE
    ) %>%
    column_spec(1, bold = TRUE) %>%
    add_header_above(c(" " = 2, "Regression Estimates" = 5, " " = 1))



# Print the table
print(publication_table)

#Save the table as an HTML file
save_kable(publication_table, 
           file = "output/tables/differences_in_treatment_groups_genes.html")


# Convert the saved HTML file to an image
webshot("output/tables/differences_in_treatment_groups_genes.html", 
        "output/tables/differences_in_treatment_groups_genes.png")

# ***********************************************************
# Title: Predicting the health outcomes of infections in hybrid mice

# Purpose: This script analysis the potential of each immune 
# gene used int the experiment to predict the 
# the weight loss during infection
# Authors: Fay Webster

field <- hm %>%
    filter(origin == "Field")

# List of dependent variables
dependent_vars <- c("IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10",
                    "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                    "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                    "TICAM1", "TNF")

# Assuming 'field' is your dataframe
#field$current_infection <- factor(field$current_infection, levels = c("uninfected", "E_falciformis", "E_ferrisi"))

# Perform regressions
results <- lapply(dependent_vars, function(var) {
    lm_formula <- as.formula(paste(var, "~ species_Eimeria", sep = " "))
    lm(lm_formula, data = field)
})

names(results) <- dependent_vars  # Name each regression with the name of the dependent variable


# Correctly add 'model' column to each tidied model's dataframe
tidy_models <- do.call(rbind, lapply(names(results), function(name) {
    model_df <- tidy(results[[name]], conf.int = TRUE)
    model_df$model <- name  # Add model name as a new column
    return(model_df)
}))

# Now 'tidy_models' should have a 'model' column correctly populated

# Filter out intercept terms
tidy_models_no_intercept <- tidy_models[!grepl("intercept", tidy_models$term, ignore.case = TRUE),]

tidy_models_no_intercept <- as.data.frame(tidy_models_no_intercept) %>%
    mutate(term = case_when(
        term == "species_EimeriaE. falciformis" ~ "E. falciformis",
        term == "species_EimeriaE. ferrisi" ~ "E. ferrisi",
    ))



# Determine the common y-axis range
common_y_limits <- range(-7,7)

# Create the coefficient plot
coef_mmr_B <- ggplot(tidy_models_no_intercept, aes(x = model, y = estimate, color = term)) +
    coord_flip() +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("E. falciformis" = "#FF0000", "E. ferrisi" = "#7A0092")) +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Coefficients estimate (Difference to uninfected)") +
    scale_y_continuous(limits = common_y_limits) +  # Apply shared y-axis limits
    theme(legend.title = element_blank(),
          legend.position = "none")


print(coef_mmr_B)

ggsave(filename = paste0(an_fi, "/coef_plot_field_genes.pdf"),
       plot = coef_mmr_B, width = 6, height = 4, dpi = 300)



####################################################
####################################################
############################# gene expression distribution
field %>%
    drop_na(species_Eimeria) %>%
    pivot_longer(cols = all_of(dependent_vars), 
                 names_to = "Genes", values_to = "Expression") %>%
    ggplot(aes(x = Expression, fill = species_Eimeria)) +
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) +
    geom_boxplot(position = "dodge2",
                 width = .5, 
                 outlier.shape = NA,
                 orientation = "y"  # Set orientation to y
    ) +
    facet_wrap(~Genes,  scales = 'free', ncol = 4) +
    labs(x = "Expression Level", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = color_mapping_f, labels = labels_f)  +
    theme(legend.title = element_blank(), 
          legend.position = c(0.85, 0.06),
          legend.text = element_markdown())+
    labs(y = "Density", 
         x = "Gene expression level") -> density_imm

density_imm

ggsave(filename = paste0(an_fi, "/density_immune_genes_field.pdf"),
       plot = density_imm, width = 10, height = 8, dpi = 300)


#######################
# combine
comb <- (density_imm | coef_mmr) +
    #  plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add fieldels (A, B, C, etc.)

# Add a figure title
comb <- comb + 
    plot_annotation(title = 'Fig. 3', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))

comb <- comb + 
    plot_layout(heights = c(1, 1), 
                widths = c(2, 1))



# Display the panel figure
# print(comb)

# Save the panel figure
ggsave(paste0(panels_fi, "/panel_immune_gene_expression_field.jpeg"), 
       comb, width = 12, height = 6, dpi = 300)  


# Update gene names to replace "." with "-" and "IFNg" with "IFNy"
regression_table <- tidy_models_no_intercept %>%
    dplyr::select(model, term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(
        model = gsub("\\.", "-", model),  # Replace "." with "-"
        model = gsub("IFNg", "IFNy", model),  # Replace "IFNg" with "IFNy"
        significance = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
        )
    ) %>%
    rename(
        Gene = model,
        'Treatment group' = term,
        Estimate = estimate,
        `Std. Error` = std.error,
        `CI Lower` = conf.low,
        `CI Upper` = conf.high,
        `P-value` = p.value,
        `Significance` = significance
    )

# Embed significance codes in the caption directly using HTML
caption_text <- "Table 2: Regression results for immune gene expression in <em>E. ferrisi</em> and <em>E. falciformis</em>-infected field-caught mice.<br><em>Significance codes:</em> '***' 0.001, '**' 0.01, '*' 0.05"

publication_table <- regression_table %>%
    kbl(
        caption = caption_text,
        format = "html",  # Or "latex" for LaTeX output
        digits = 3,
        escape = FALSE  # allow HTML in the caption
    ) %>%
    kable_styling(
        bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
        full_width = FALSE
    ) %>%
    column_spec(1, bold = TRUE) %>%
    add_header_above(c(" " = 2, "Regression Estimates" = 5, " " = 1))


# Print the table
print(publication_table)

#Save the table as an HTML file
save_kable(publication_table, 
           file = "output/tables/differences_in_treatment_groups_genes_field.html")


# Convert the saved HTML file to an image
webshot("output/tables/differences_in_treatment_groups_genes_field.html", 
        "output/tables/differences_in_treatment_groups_genes_field.png")

# Create simplified table with only significant results
create_simplified_table <- function(tidy_data, dataset_name) {
    significant_results <- tidy_data %>%
        dplyr::filter(p.value < 0.05) %>%
        dplyr::select(model, term, estimate, p.value) %>%
        mutate(
            estimate = round(estimate, 2),
            p.value = case_when(
                p.value < 0.001 ~ "< 0.001",
                p.value < 0.01 ~ paste("=", round(p.value, 3)),
                TRUE ~ paste("=", round(p.value, 2))
            )
        ) %>%
        rename(Gene = model, Species = term, Estimate = estimate, `P-value` = p.value) %>%
        arrange(Gene, Species)
    
    return(significant_results)
}

# Create simplified tables
lab_simple <- create_simplified_table(lab_analysis$tidy_results, "Laboratory")
field_simple <- create_simplified_table(field_analysis$tidy_results, "Field")

lab_simple
field_simple

# Create publication-ready simplified tables
create_publication_simple_table <- function(simple_data, caption_text, filename) {
    
    # Format the table
    pub_table <- simple_data %>%
        kbl(caption = caption_text, 
            format = "html", 
            digits = 2,
            escape = FALSE) %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                      full_width = FALSE,
                      font_size = 12) %>%
        column_spec(1, bold = TRUE) %>%
        column_spec(2, italic = TRUE) %>%  # Italicize species names
        column_spec(3, width = "2cm") %>%
        column_spec(4, width = "2cm") %>%
        add_header_above(c(" " = 2, "Regression Results" = 2))
    
    # Save table
    save_kable(pub_table, file = filename)
    
    return(pub_table)
}

# Create the publication tables
lab_pub_simple <- create_publication_simple_table(
    lab_simple,
    "Table 1: Significant immune gene expression changes in laboratory-infected mice",
    "output/tables/lab_genes_simple.html"
)

field_pub_simple <- create_publication_simple_table(
    field_simple, 
    "Table 2: Significant immune gene expression changes in field-caught mice",
    "output/tables/field_genes_simple.html"
)

# Print the tables
print(lab_pub_simple)
print(field_pub_simple)

# Optional: Create a combined table showing both lab and field side by side
combined_simple <- lab_simple %>%
    rename(Lab_Estimate = Estimate, Lab_Pvalue = `P-value`) %>%
    full_join(field_simple %>% 
                  rename(Field_Estimate = Estimate, Field_Pvalue = `P-value`),
              by = c("Gene", "Species")) %>%
    arrange(Gene, Species)

# Convert to GT table
combined_gt_table <- combined_simple %>%
    gt() %>%
    # Add title and subtitle
    tab_header(
        title = "Comparison of significant immune gene responses",
        subtitle = "Controlled laboratory infections vs. field-caught mice infected with Eimeria spp."
    ) %>%
    # Format column headers
    cols_label(
        Gene = "Gene",
        Species = "Species", 
        Lab_Estimate = "Estimate",
        Lab_Pvalue = "P-value",
        Field_Estimate = "Estimate", 
        Field_Pvalue = "P-value"
    ) %>%
    # Add column spanners
    tab_spanner(
        label = "Controlled laboratory infections",
        columns = c(Lab_Estimate, Lab_Pvalue)
    ) %>%
    tab_spanner(
        label = "Field-caught mice",
        columns = c(Field_Estimate, Field_Pvalue)
    ) %>%
    # Style the table
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = Gene)
    ) %>%
    tab_style(
        style = cell_text(style = "italic"),
        locations = cells_body(columns = Species)
    ) %>%
    # Format numbers
    fmt_number(
        columns = c(Lab_Estimate, Field_Estimate),
        decimals = 2
    ) %>%
    # Handle missing values
    sub_missing(
        columns = everything(),
        missing_text = "—"
    ) %>%
    # Add source note with italicized Eimeria
    tab_source_note(
        source_note = md("Only genes with p < 0.05 in at least one condition shown. Full results in Supplementary Tables S1-S2. Infections with *Eimeria* spp.")
    ) %>%
    # Style the table
    tab_options(
        table.font.size = 12,
        data_row.padding = px(3),
        summary_row.padding = px(3),
        grand_summary_row.padding = px(3),
        footnotes.padding = px(3),
        source_notes.padding = px(3),
        row_group.padding = px(3)
    )

# If you want to also italicize Eimeria in the subtitle:
combined_gt_table <- combined_gt_table %>%
    tab_header(
        title = "Comparison of significant immune gene responses",
        subtitle = md("Controlled laboratory infections vs. field-caught mice infected with *Eimeria* spp.")
    )

# Print the table
print(combined_gt_table)

# Save using your function
save_table_all_formats(combined_gt_table, "gene_responses_combined_gt", "output/tables")
