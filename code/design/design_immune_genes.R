# Create Supplementary Table for Immune Genes
# Multiple approaches provided


# Create the data frame
immune_genes_table <- data.frame(
    Gene_Symbol = c(
        # Mucosal Defense
        "RETNLB", "MUC2", "MUC5AC",
        # Innate Immune Markers  
        "NCR1", "MPO",
        # Inflammatory Signaling
        "MYD88", "CASP1", "IL1RN", "TNF", "IL.6",
        # Chemokines and Migration
        "CXCL9", "CXCR3", 
        # Th1/Th2 Response
        "IFNy", "IL.13", "IRGM1",
        # Regulatory Pathways
        "TICAM1", "SOCS1", "IDO1",
        # Cytotoxic Markers
        "PRF1"
    ),
    
    Full_Gene_Name = c(
        # Mucosal Defense
        "Resistin-like beta", "Mucin 2", "Mucin 5AC",
        # Innate Immune Markers
        "Natural killer cell p46-related protein", "Myeloperoxidase",
        # Inflammatory Signaling  
        "Myeloid differentiation primary response 88", "Caspase 1", 
        "Interleukin 1 receptor antagonist", "Tumor necrosis factor alpha", "Interleukin 6",
        # Chemokines and Migration
        "C-X-C motif chemokine ligand 9", "C-X-C motif chemokine receptor 3",
        # Th1/Th2 Response
        "Interferon gamma", "Interleukin 13", "Immunity-related GTPase M1",
        # Regulatory Pathways
        "TIR domain-containing adaptor molecule 1", "Suppressor of cytokine signaling 1", 
        "Indoleamine 2,3-dioxygenase 1",
        # Cytotoxic Markers
        "Perforin 1"
    ),
    
    Functional_Category = c(
        # Mucosal Defense
        "Mucosal barrier", "Mucosal barrier", "Mucosal barrier",
        # Innate Immune Markers
        "NK cell marker", "Neutrophil marker", 
        # Inflammatory Signaling
        "TLR signaling", "Inflammasome", "Anti-inflammatory", 
        "Pro-inflammatory cytokine", "Pro-inflammatory cytokine",
        # Chemokines and Migration
        "Chemokine", "Chemokine receptor",
        # Th1/Th2 Response  
        "Th1 effector cytokine", "Th2 cytokine", "Intracellular defense",
        # Regulatory Pathways
        "TLR signaling", "Immune regulation", "Immune regulation",
        # Cytotoxic Markers
        "Cytotoxic effector"
    ),
    
    Primary_Role = c(
        # Mucosal Defense
        "Goblet cell defense factor, mucus production",
        "Major secretory mucin in gastrointestinal tract", 
        "Surface goblet cell mucin production",
        # Innate Immune Markers
        "Natural killer cell activation and cytotoxicity",
        "Neutrophil-mediated oxidative pathogen killing",
        # Inflammatory Signaling
        "TLR-mediated NF-κB activation, inflammation",
        "IL-1β and IL-18 production, pyroptosis",
        "Natural IL-1β antagonist, infection control",
        "Inflammation, immune cell activation", 
        "Inflammation regulation, acute phase response",
        # Chemokines and Migration
        "Immune cell migration, Th1 activation",
        "CXCL9/CXCL11 receptor, T cell trafficking",
        # Th1/Th2 Response
        "Macrophage activation, pathogen clearance",
        "Alternative macrophage activation, mucus production",
        "Autonomous cell defense, parasite clearance",
        # Regulatory Pathways
        "Type I IFN production, antiviral response",
        "JAK/STAT pathway regulation, T-cell differentiation",
        "T-cell activity regulation via tryptophan depletion",
        # Cytotoxic Markers
        "Cytotoxic granule protein, target cell lysis"
    ),
    
    Category_Group = c(
        # Mucosal Defense
        rep("Mucosal Defense", 3),
        # Innate Immune Markers
        rep("Innate Immune Markers", 2),
        # Inflammatory Signaling
        rep("Inflammatory Signaling", 5),
        # Chemokines and Migration
        rep("Chemokines and Migration", 2),
        # Th1/Th2 Response
        rep("Th1/Th2 Response", 3),
        # Regulatory Pathways
        rep("Regulatory Pathways", 3),
        # Cytotoxic Markers
        rep("Cytotoxic Markers", 1)
    )
)

# Method 1: Basic kable table
basic_table <- kable(immune_genes_table[,1:4], 
                     caption = "Supplementary Table X: Immune gene panel used for expression analysis",
                     col.names = c("Gene Symbol", "Full Gene Name", "Functional Category", "Primary Role"))

print(basic_table)

# Method 2: Enhanced kableExtra table with grouping
enhanced_table <- immune_genes_table %>%
    dplyr::select(-Category_Group) %>%
    kable(caption = "Immune gene panel used for expression analysis",
          col.names = c("Gene Symbol", "Full Gene Name", "Functional Category", "Primary Role")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    pack_rows("Mucosal Defense", 1, 3) %>%
    pack_rows("Innate Immune Markers", 4, 5) %>%
    pack_rows("Inflammatory Signaling", 6, 10) %>%
    pack_rows("Chemokines and Migration", 11, 12) %>%
    pack_rows("Th1/Th2 Response", 13, 15) %>%
    pack_rows("Regulatory Pathways", 16, 18) %>%
    pack_rows("Cytotoxic Markers", 19, 19)

enhanced_table



# Method 3: Save to your output directory structure
write.csv(immune_genes_table, "output/tables/Supplementary_table_immune_genes.csv", row.names = FALSE)

# Also save the basic kable version
writeLines(as.character(basic_table), "output/tables/Supplementary_table_immune_genes.html")

# Save enhanced version if it works
# enhanced_table %>% save_kable("output/tables/Supplementary_table_immune_genes_enhanced.html")

# Method 4: Create Word document table (requires flextable)
# install.packages("flextable") # if needed
library(flextable)

word_table <- immune_genes_table %>%
    dplyr::select(-Category_Group) %>%
    flextable() %>%
    set_header_labels(
        Gene_Symbol = "Gene Symbol",
        Full_Gene_Name = "Full Gene Name", 
        Functional_Category = "Functional Category",
        Primary_Role = "Primary Role"
    ) %>%
    theme_vanilla() %>%
    autofit() %>%
    add_header_row(
        values = c("Supplementary Table X: Immune gene panel used for expression analysis"),
        colwidths = 4
    )

# Save as Word document
save_as_docx(word_table, path = "output/tables/Supplementary_table_immune_genes.docx")

# Method 5: For LaTeX/PDF output
latex_table <- immune_genes_table %>%
    dplyr::select(-Category_Group) %>%
    kable(format = "latex", 
          caption = "Immune gene panel used for expression analysis",
          col.names = c("Gene Symbol", "Full Gene Name", "Functional Category", "Primary Role"),
          booktabs = TRUE) %>%
    kable_styling(latex_options = c("striped", "scale_down"))

# Save LaTeX table
writeLines(as.character(latex_table), "output/tables/Supplementary_table_immune_genes.tex")

# Print the genes in order for verification
cat("Genes in order:\n")
cat(paste(immune_genes_table$Gene_Symbol, collapse = ", "))

# Verify against your gene list
your_genes <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", 
                "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

# Check if all genes are included
missing_genes <- setdiff(your_genes, immune_genes_table$Gene_Symbol)
extra_genes <- setdiff(immune_genes_table$Gene_Symbol, your_genes)

cat("\nMissing genes:", missing_genes)
cat("\nExtra genes:", extra_genes)

# Create Supplementary Table for Immune Genes
# Multiple approaches provided

library(knitr)
library(kableExtra)
library(dplyr)

# Create the data frame
immune_genes_table <- data.frame(
    Gene_Symbol = c(
        # Mucosal Defense
        "RETNLB", "MUC2", "MUC5AC",
        # Innate Immune Markers  
        "NCR1", "MPO",
        # Inflammatory Signaling
        "MYD88", "CASP1", "IL1RN", "TNF", "IL.6",
        # Chemokines and Migration
        "CXCL9", "CXCR3", 
        # Th1/Th2 Response
        "IFNy", "IL.13", "IRGM1",
        # Regulatory Pathways
        "TICAM1", "SOCS1", "IDO1",
        # Cytotoxic Markers
        "PRF1"
    ),
    
    Full_Gene_Name = c(
        # Mucosal Defense
        "Resistin-like beta", "Mucin 2", "Mucin 5AC",
        # Innate Immune Markers
        "Natural killer cell p46-related protein", "Myeloperoxidase",
        # Inflammatory Signaling  
        "Myeloid differentiation primary response 88", "Caspase 1", 
        "Interleukin 1 receptor antagonist", "Tumor necrosis factor alpha", "Interleukin 6",
        # Chemokines and Migration
        "C-X-C motif chemokine ligand 9", "C-X-C motif chemokine receptor 3",
        # Th1/Th2 Response
        "Interferon gamma", "Interleukin 13", "Immunity-related GTPase M1",
        # Regulatory Pathways
        "TIR domain-containing adaptor molecule 1", "Suppressor of cytokine signaling 1", 
        "Indoleamine 2,3-dioxygenase 1",
        # Cytotoxic Markers
        "Perforin 1"
    ),
    
    Functional_Category = c(
        # Mucosal Defense
        "Mucosal barrier", "Mucosal barrier", "Mucosal barrier",
        # Innate Immune Markers
        "NK cell marker", "Neutrophil marker", 
        # Inflammatory Signaling
        "TLR signaling", "Inflammasome", "Anti-inflammatory", 
        "Pro-inflammatory cytokine", "Pro-inflammatory cytokine",
        # Chemokines and Migration
        "Chemokine", "Chemokine receptor",
        # Th1/Th2 Response  
        "Th1 effector cytokine", "Th2 cytokine", "Intracellular defense",
        # Regulatory Pathways
        "TLR signaling", "Immune regulation", "Immune regulation",
        # Cytotoxic Markers
        "Cytotoxic effector"
    ),
    
    Primary_Role = c(
        # Mucosal Defense
        "Goblet cell defense factor, mucus production",
        "Major secretory mucin in gastrointestinal tract", 
        "Surface goblet cell mucin production",
        # Innate Immune Markers
        "Natural killer cell activation and cytotoxicity",
        "Neutrophil-mediated oxidative pathogen killing",
        # Inflammatory Signaling
        "TLR-mediated NF-κB activation, inflammation",
        "IL-1β and IL-18 production, pyroptosis",
        "Natural IL-1β antagonist, infection control",
        "Inflammation, immune cell activation", 
        "Inflammation regulation, acute phase response",
        # Chemokines and Migration
        "Immune cell migration, Th1 activation",
        "CXCL9/CXCL11 receptor, T cell trafficking",
        # Th1/Th2 Response
        "Macrophage activation, pathogen clearance",
        "Alternative macrophage activation, mucus production",
        "Autonomous cell defense, parasite clearance",
        # Regulatory Pathways
        "Type I IFN production, antiviral response",
        "JAK/STAT pathway regulation, T-cell differentiation",
        "T-cell activity regulation via tryptophan depletion",
        # Cytotoxic Markers
        "Cytotoxic granule protein, target cell lysis"
    ),
    
    Category_Group = c(
        # Mucosal Defense
        rep("Mucosal Defense", 3),
        # Innate Immune Markers
        rep("Innate Immune Markers", 2),
        # Inflammatory Signaling
        rep("Inflammatory Signaling", 5),
        # Chemokines and Migration
        rep("Chemokines and Migration", 2),
        # Th1/Th2 Response
        rep("Th1/Th2 Response", 3),
        # Regulatory Pathways
        rep("Regulatory Pathways", 3),
        # Cytotoxic Markers
        rep("Cytotoxic Markers", 1)
    )
)

# Method 1: Basic kable table
basic_table <- kable(immune_genes_table[,1:4], 
                     caption = "Supplementary Table X: Immune gene panel used for expression analysis",
                     col.names = c("Gene Symbol", "Full Gene Name", "Functional Category", "Primary Role"))

print(basic_table)

# Method 2: Enhanced kableExtra table with grouping
enhanced_table <- immune_genes_table %>%
    dplyr::select(-Category_Group) %>%
    kable(caption = "Supplementary Table X: Immune gene panel used for expression analysis",
          col.names = c("Gene Symbol", "Full Gene Name", "Functional Category", "Primary Role")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    pack_rows("Mucosal Defense", 1, 3) %>%
    pack_rows("Innate Immune Markers", 4, 5) %>%
    pack_rows("Inflammatory Signaling", 6, 10) %>%
    pack_rows("Chemokines and Migration", 11, 12) %>%
    pack_rows("Th1/Th2 Response", 13, 15) %>%
    pack_rows("Regulatory Pathways", 16, 18) %>%
    pack_rows("Cytotoxic Markers", 19, 19)

enhanced_table

enhanced_table

# Method 3: Save to your output directory structure
write.csv(immune_genes_table, "output/tables/Supplementary_table_immune_genes.csv", row.names = FALSE)

# Convert table to a ggplot object for better image control
table_gg <- tableGrob(immune_genes_table[,1:4], rows = NULL)

# Add title
title <- textGrob("Supplementary Table X: Immune gene panel used for expression analysis", 
                  gp = gpar(fontsize = 14, fontface = "bold"))

# Combine title and table
combined <- arrangeGrob(title, table_gg, heights = c(0.1, 0.9))

# Save as JPEG with high quality
ggsave("output/tables/Supplementary_table_immune_genes_ggplot.jpeg", 
       combined, 
       width = 12, height = 10, 
       dpi = 300, 
       device = "jpeg")

# Also save as PNG for better quality
ggsave("output/tables/Supplementary_table_immune_genes_ggplot.png", 
       combined, 
       width = 12, height = 10, 
       dpi = 300, 
       device = "png")


