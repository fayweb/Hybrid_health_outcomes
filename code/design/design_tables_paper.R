# Function to save gt tables in multiple formats
save_table_all_formats <- function(table_object, table_name, output_dir = tables) {
    
    # Create base filename
    base_path <- file.path(output_dir, table_name)
    
    # Save in all formats
    tryCatch({
        # HTML (best for Google Docs)
        gtsave(table_object, paste0(base_path, ".html"))
        cat("✓ Saved", table_name, "as HTML\n")
        
        # Word document
        gtsave(table_object, paste0(base_path, ".docx"))
        cat("✓ Saved", table_name, "as DOCX\n")
        
        # PNG (high quality image)
        gtsave(table_object, paste0(base_path, ".png"), 
               vwidth = 1000, vheight = 800)
        cat("✓ Saved", table_name, "as PNG\n")
        
        # PDF (best for LaTeX)
        gtsave(table_object, paste0(base_path, ".pdf"))
        cat("✓ Saved", table_name, "as PDF\n")
        
        # LaTeX code
        table_object %>%
            as_latex() %>%
            writeLines(paste0(base_path, ".tex"))
        cat("✓ Saved", table_name, "as TEX\n")
        
        cat("✅ All formats saved successfully for", table_name, "\n\n")
        
    }, error = function(e) {
        cat("❌ Error saving", table_name, ":", e$message, "\n")
    })
}




# =============================================
# TABLE 1: STRAIN COMPOSITION AND ORIGINS
# =============================================

# Step 1: UNGROUP the data first
lab_clean <- lab %>% ungroup()

# Step 2: Simple strain count
strain_counts <- lab_clean %>%
    count(mouse_strain) %>%
    arrange(desc(n))

print(strain_counts)

# Step 3: Create a clean, organized table
table1_data <- lab_clean %>%
    mutate(
        strain_type = case_when(
            mouse_strain == "SCHUNT SCHUNT" ~ "Pure M. m. domesticus",
            mouse_strain == "PWD PWD" ~ "Pure M. m. musculus", 
            mouse_strain == "STRA STRA" ~ "Pure M. m. domesticus",
            mouse_strain == "BUSNA BUSNA" ~ "Pure M. m. musculus",
            mouse_strain == "NMRI" ~ "Laboratory strain",
            TRUE ~ "F1 hybrid"
        ),
        strain_name = case_when(
            mouse_strain == "SCHUNT SCHUNT" ~ "SCHUNT",
            mouse_strain == "PWD PWD" ~ "PWD",
            mouse_strain == "STRA STRA" ~ "STRA", 
            mouse_strain == "BUSNA BUSNA" ~ "BUSNA",
            mouse_strain == "NMRI" ~ "NMRI",
            TRUE ~ as.character(mouse_strain)
        )
    ) %>%
    count(strain_type, strain_name) %>%
    arrange(strain_type, desc(n))

print(table1_data)

# Step 4: Create the final publication-ready table
table1 <- table1_data %>%
    gt() %>%
    tab_header(
        title = "Table 1. Mouse strain composition and genetic background",
        subtitle = "Distribution of 136 mice across subspecies and strain types"
    ) %>%
    cols_label(
        strain_type = "Subspecies/Type",
        strain_name = "Strain",
        n = "n"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    tab_footnote(
        footnote = "Origins: SCHUNT (Schweben, Germany), PWD (Kunratice, Czech Republic), STRA (Straas, Germany), BUSNA (Buškovice, Czech Republic)",
        locations = cells_column_labels(columns = strain_type)
    ) %>%
    # Add subtotals by strain type
    summary_rows(
        groups = TRUE,
        columns = n,
        fns = list(Subtotal = ~sum(.))
    )

print(table1)

save_table_all_formats(table1, "Table1_strain_composition")



# =============================================
# TABLE 2: EXPERIMENTAL DESIGN OVERVIEW  
# =============================================

# =============================================
# TABLE 2: EXPERIMENTAL DESIGN OVERVIEW  
# =============================================

table2 <- lab_clean %>%
    count(experiment, current_infection) %>%
    arrange(experiment, current_infection) %>%
    pivot_wider(
        names_from = current_infection, 
        values_from = n, 
        values_fill = 0
    ) %>%
    # Add total column
    rowwise() %>%
    mutate(Total = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    ungroup() %>%
    # Add total row manually
    bind_rows(
        summarise(., 
                  experiment = "TOTAL",
                  across(where(is.numeric), sum)
        )
    ) %>%
    gt() %>%
    tab_header(
        title = "Table 2. Experimental design and infection protocols",
        subtitle = "Sample distribution across experiments and final infection status"
    ) %>%
    cols_label(
        experiment = "Experiment"
    ) %>%
    tab_spanner(
        label = "Final Infection Status (n mice)",
        columns = c("E. ferrisi", "E. falciformis", "Uninfected controls")
    ) %>%
    # Make species names italic in column headers
    cols_label(
        "E. ferrisi" = html("<em>E. ferrisi</em>"),
        "E. falciformis" = html("<em>E. falciformis</em>")
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Style the total row
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(rows = experiment == "TOTAL")
    )

print(table2)

save_tables <- function(x) {
    
}


save_table_all_formats(table2, "Table2_experimental_design")


# =============================================
# TABLE 3: INFECTION HISTORY BREAKDOWN
# =============================================

infection_breakdown <- lab_clean %>%
    count(infection_history, immunization, current_infection) %>%
    arrange(infection_history)

table3 <- infection_breakdown %>%
    mutate(
        primary_infection = str_extract(infection_history, "^[^_]+"),
        challenge_infection = str_extract(infection_history, "[^_]+$"),
        primary_infection = case_when(
            primary_infection == "ferrisi" ~ "E. ferrisi",
            primary_infection == "falciformis" ~ "E. falciformis", 
            primary_infection == "uninfected" ~ "Uninfected",
            TRUE ~ primary_infection
        ),
        challenge_infection = case_when(
            challenge_infection == "ferrisi" ~ "E. ferrisi",
            challenge_infection == "falciformis" ~ "E. falciformis",
            challenge_infection == "uninfected" ~ "Uninfected", 
            TRUE ~ challenge_infection
        ),
        # Add HTML formatting for italics in the data
        primary_infection = case_when(
            str_detect(primary_infection, "E\\.") ~ paste0("<em>", primary_infection, "</em>"),
            TRUE ~ primary_infection
        ),
        challenge_infection = case_when(
            str_detect(challenge_infection, "E\\.") ~ paste0("<em>", challenge_infection, "</em>"),
            TRUE ~ challenge_infection
        ),
        current_infection = case_when(
            str_detect(current_infection, "E\\.") ~ paste0("<em>", str_replace(current_infection, "E\\. ", "E. "), "</em>"),
            TRUE ~ as.character(current_infection)
        )
    ) %>%
    .[, c("primary_infection", "challenge_infection", "immunization", "current_infection", "n")] %>%
    gt() %>%
    tab_header(
        title = "Table 3. Infection history and immunization status",
        subtitle = "Primary and challenge infection combinations"
    ) %>%
    cols_label(
        primary_infection = "Primary Infection",
        challenge_infection = "Challenge Infection", 
        immunization = "Immune Status",
        current_infection = "Final Status",
        n = "n"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Enable HTML rendering
    fmt_markdown(columns = c(primary_infection, challenge_infection, current_infection))

print(table3)

save_table_all_formats(table3, "Table3_infection_history")
# =============================================
# TABLE 4: STRAIN × INFECTION BREAKDOWN
# =============================================

strain_infection <- lab_clean %>%
    mutate(
        strain_group = case_when(
            mouse_strain %in% c("SCHUNT SCHUNT", "STRA STRA") ~ "M. m. domesticus",
            mouse_strain %in% c("PWD PWD", "BUSNA BUSNA") ~ "M. m. musculus",
            mouse_strain == "NMRI" ~ "NMRI",
            TRUE ~ "F1 hybrids"
        )
    ) %>%
    count(strain_group, current_infection) %>%
    pivot_wider(names_from = current_infection, values_from = n, values_fill = 0) %>%
    mutate(Total = rowSums(across(where(is.numeric)))) %>%
    # Add total row manually to avoid summary_rows error
    bind_rows(
        summarise(., 
                  strain_group = "TOTAL",
                  across(where(is.numeric), sum)
        )
    )

table4 <- strain_infection %>%
    gt() %>%
    tab_header(
        title = "Strain distribution across infection groups",
        subtitle = "Final infection status by genetic background"
    ) %>%
    cols_label(
        strain_group = "Strain Group",
        "E. ferrisi" = html("<em>E. ferrisi</em>"),
        "E. falciformis" = html("<em>E. falciformis</em>")
    ) %>%
    tab_spanner(
        label = "Final Infection Status",
        columns = c("E. ferrisi", "E. falciformis", "Uninfected controls")
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Style the total row
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(rows = strain_group == "TOTAL")
    )

print(table4)


save_table_all_formats(table4, "Supplementary_table_strain_infection_groups")

# =============================================
# SUPPLEMENTARY TABLE: DETAILED CROSSINGS
# =============================================

detailed_crossings <- lab_clean %>%
    filter(mouse_strain != "NMRI") %>%
    count(mouse_strain, hybrid_status, current_infection) %>%
    arrange(mouse_strain, current_infection)

supp_table <- detailed_crossings %>%
    pivot_wider(names_from = current_infection, values_from = n, values_fill = 0) %>%
    mutate(
        Total = rowSums(across(where(is.numeric))),
        # Add HTML formatting for italics in the hybrid_status column
        hybrid_status = case_when(
            str_detect(hybrid_status, "M\\. m\\.") ~ str_replace_all(hybrid_status, 
                                                                     "(M\\. m\\. \\w+)", "<em>\\1</em>"),
            TRUE ~ hybrid_status
        )
    ) %>%
    gt() %>%
    tab_header(
        title = "Detailed strain crosses and infections",
        subtitle = "Complete breakdown of all wild-derived strains and hybrids"
    ) %>%
    cols_label(
        mouse_strain = "Strain Cross",
        hybrid_status = "Genetic Background",
        "E. ferrisi" = html("<em>E. ferrisi</em>"),
        "E. falciformis" = html("<em>E. falciformis</em>")
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    # Enable HTML rendering for the hybrid_status column
    fmt_markdown(columns = hybrid_status)

print(supp_table)
save_table_all_formats(supp_table, "Supplementary_table_detailed_crossings")

