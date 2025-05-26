# =============================================================================
# COMPLETE WORKING FIELD DATA VISUALIZATIONS
# =============================================================================
# All your data is perfect! Let's create the full visualization suite.

library(ggplot2)
library(dplyr)
library(patchwork)

# =============================================================================
# 1. DETAILED SAMPLE FLOWCHART
# =============================================================================

create_detailed_flowchart <- function(field) {
    flowchart_data <- data.frame(
        step = factor(c("Total Field\nMice", "With Immune\nData", "Known Infection\nStatus", 
                        "Species\nClassified", "Intensity\nData", "Complete\nData"),
                      levels = c("Total Field\nMice", "With Immune\nData", "Known Infection\nStatus", 
                                 "Species\nClassified", "Intensity\nData", "Complete\nData")),
        n = c(336, 336, 305, 169, 92, 14)  # Your exact numbers
    )
    
    p <- ggplot(flowchart_data, aes(x = step, y = n)) +
        geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
        geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
        labs(title = "Sample Size Progression Through Analysis Pipeline",
             subtitle = "From field collection to analysis-ready datasets",
             x = "Analysis Stage", y = "Number of Mice") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            plot.title = element_text(size = 14, face = "bold"),
            panel.grid.major.x = element_blank()
        ) +
        ylim(0, 350)
    
    return(p)
}

# =============================================================================
# 2. METHOD AVAILABILITY BAR CHART
# =============================================================================

create_method_availability <- function(field) {
    method_data <- data.frame(
        Method = c("Immune Genes", "Tissue qPCR", "Feces qPCR", "Oocyst Count", "Amplicon Seq", "Species ID"),
        Count = c(336, 185, 156, 179, 134, 169),  # Your exact numbers
        Percentage = round(c(336, 185, 156, 179, 134, 169) / 336 * 100, 1)
    ) %>%
        mutate(Method = factor(Method, levels = rev(Method)))  # Reverse for nice display
    
    p <- ggplot(method_data, aes(x = Method, y = Count)) +
        geom_col(fill = "forestgreen", alpha = 0.8) +
        geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
                  hjust = -0.1, size = 3.5) +
        coord_flip() +
        labs(title = "Data Availability by Quantification Method",
             subtitle = "Number and percentage of mice with each data type",
             x = "Quantification Method", y = "Number of Mice") +
        theme_minimal() +
        theme(plot.title = element_text(size = 14, face = "bold")) +
        ylim(0, 360)
    
    return(p)
}

# =============================================================================
# 3. INFECTION STATUS PIE CHART
# =============================================================================

create_infection_pie <- function(field) {
    status_data <- data.frame(
        status = c("Infected", "Uninfected"),
        n = c(133, 172),  # Your exact numbers
        percentage = c(43.6, 56.4)
    ) %>%
        mutate(label = paste0(status, "\n", n, " mice\n(", percentage, "%)"))
    
    p <- ggplot(status_data, aes(x = "", y = n, fill = status)) +
        geom_col(color = "white", size = 1) +
        coord_polar(theta = "y") +
        geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
        labs(title = "Infection Status Distribution",
             subtitle = "Based on 305 mice with known status") +
        scale_fill_manual(values = c("Infected" = "#d73027", "Uninfected" = "#4575b4")) +
        theme_void() +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
    
    return(p)
}

# =============================================================================
# 4. SPECIES DISTRIBUTION
# =============================================================================

create_species_bar <- function(field) {
    species_data <- data.frame(
        species = c("Uninfected", "E. ferrisi", "E. falciformis"),
        n = c(110, 45, 14),  # Your exact numbers
        percentage = round(c(110, 45, 14) / 169 * 100, 1)
    ) %>%
        mutate(species = factor(species, levels = species))
    
    p <- ggplot(species_data, aes(x = species, y = n, fill = species)) +
        geom_col(alpha = 0.8, width = 0.7) +
        geom_text(aes(label = paste0(n, "\n(", percentage, "%)")), 
                  vjust = -0.5, size = 4) +
        labs(title = "Eimeria Species Distribution",
             subtitle = "Based on 169 mice with species classification",
             x = "Species", y = "Number of Mice") +
        scale_fill_manual(values = c("Uninfected" = "#4575b4", "E. ferrisi" = "#fee08b", 
                                     "E. falciformis" = "#d73027")) +
        theme_minimal() +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 14, face = "bold"),
            panel.grid.major.x = element_blank()
        ) +
        ylim(0, 120)
    
    return(p)
}

# =============================================================================
# 5. INFECTION INTENSITY BY SPECIES
# =============================================================================

create_intensity_boxplot <- function(field) {
    # Filter for mice with both species ID and qPCR data
    intensity_data <- field %>%
        filter(has_tissue_qPCR & species_clean %in% c("E. ferrisi", "E. falciformis", "Uninfected")) %>%
        filter(!is.na(delta_ct_cewe_MminusE))
    
    p <- ggplot(intensity_data, aes(x = species_clean, y = delta_ct_cewe_MminusE, fill = species_clean)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
        labs(title = "Infection Intensity by Species (ΔCt)",
             subtitle = "Lower ΔCt = Higher infection intensity",
             x = "Species", y = "ΔCt (Mouse - Eimeria)") +
        scale_fill_manual(values = c("Uninfected" = "#4575b4", "E. ferrisi" = "#fee08b", 
                                     "E. falciformis" = "#d73027")) +
        theme_minimal() +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    return(p)
}

# =============================================================================
# 6. METHOD COMBINATIONS
# =============================================================================

create_method_combos <- function(field) {
    # Your top method combinations
    combo_data <- data.frame(
        combo = c("TFI", "OASI", "TFOSI", "TI", "OI", "TOSI", "TFOASI", "FI"),
        count = c(115, 113, 22, 21, 17, 13, 10, 9),
        description = c("Tissue + Feces qPCR\n+ Immune", 
                        "Oocyst + Amplicon + Species\n+ Immune",
                        "All methods\nexcept Amplicon",
                        "Tissue qPCR\n+ Immune only",
                        "Oocyst\n+ Immune only",
                        "Tissue + Oocyst + Species\n+ Immune",
                        "All methods\ncomplete",
                        "Feces qPCR\n+ Immune only")
    ) %>%
        mutate(
            percentage = round(count / 336 * 100, 1),
            description = factor(description, levels = rev(description))
        )
    
    p <- ggplot(combo_data, aes(x = description, y = count)) +
        geom_col(fill = "darkorange", alpha = 0.8) +
        geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
                  hjust = -0.1, size = 3) +
        coord_flip() +
        labs(title = "Top Method Combinations",
             subtitle = "T=Tissue qPCR, F=Feces qPCR, O=Oocyst, A=Amplicon, S=Species, I=Immune",
             x = "Method Combination", y = "Number of Mice") +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, face = "bold"))
    
    return(p)
}

# =============================================================================
# 7. CREATE ALL PLOTS
# =============================================================================

# Create individual plots
p1_detailed <- create_detailed_flowchart(field)
p2_methods <- create_method_availability(field)
p3_infection <- create_infection_pie(field)
p4_species <- create_species_bar(field)
p5_intensity <- create_intensity_boxplot(field)
p6_combos <- create_method_combos(field)

# Create a comprehensive overview
overview_layout <- "
AABBCC
DDEEFF
"

field_overview <- p1_detailed + p2_methods + p3_infection + 
    p4_species + p5_intensity + p6_combos +
    plot_layout(design = overview_layout) +
    plot_annotation(
        title = "Field Data Analysis Overview: Wild Mouse Eimeria Infections",
        subtitle = "336 mice with immune gene expression data from Brandenburg, Germany (2016-2019)\nOverall infection rate: 43.6% based on multiple detection methods",
        caption = "Data integrated from qPCR, amplicon sequencing, and oocyst counting",
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

# =============================================================================
# 8. DISPLAY AND SAVE
# =============================================================================

cat("=== PLOTS CREATED SUCCESSFULLY ===\n")
cat("Available plots:\n")
cat("- p1_detailed: Sample size flowchart\n")
cat("- p2_methods: Method availability\n")
cat("- p3_infection: Infection status pie chart\n") 
cat("- p4_species: Species distribution\n")
cat("- p5_intensity: Infection intensity by species\n")
cat("- p6_combos: Method combinations\n")
cat("- field_overview: Complete overview (all plots combined)\n\n")

cat("To view plots:\n")
cat("field_overview  # Complete overview\n")
cat("p1_detailed     # Individual plots\n\n")

cat("To save plots:\n")
cat("ggsave('field_overview_complete.png', field_overview, width = 18, height = 12, dpi = 300)\n")
cat("ggsave('sample_flowchart.png', p1_detailed, width = 12, height = 8, dpi = 300)\n")
cat("ggsave('infection_status.png', p3_infection, width = 8, height = 8, dpi = 300)\n")

# Summary statistics for your thesis
cat("\n=== SUMMARY FOR THESIS ===\n")
cat("Field validation dataset: 336 mice with immune gene expression\n")
cat("Infection rate: 43.6% (133/305 mice with known status)\n")
cat("Species breakdown: E. ferrisi (45), E. falciformis (14), Uninfected (110)\n")
cat("Quantification methods: Tissue qPCR (185), Amplicon (134), Oocyst (179)\n")
cat("Analysis-ready subsets: Intensity data (92), Complete data (14)\n")

