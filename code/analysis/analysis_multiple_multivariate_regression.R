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



# Create the coefficient plot
ggplot(tidy_models_no_intercept, aes(x = model, y = estimate, color = term)) +
    coord_flip() +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("E. falciformis" = "salmon", "E. ferrisi" = "forestgreen")) +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Coefficients estimate (Difference to uninfected)") +
    theme(legend.title = element_blank(),
          legend.position = "none") -> coef_mmr

print(coef_mmr)

ggsave(filename = paste0(an_fi, "/coef_plot_lab_genes.jpeg"),
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
    
  rm(coef_mmr, comb, density_imm, results, tidy_models, 
     tidy_models_no_intercept, biplot, coefs5, contr_PC1, contr_PC2, 
     figure_panel, pc1_current_infection, pc2_current_infection,
     pc1_WL_current_infection, pc2_WL_current_infection, pca_individuals, vpg,
     coefs6, model_6, panel_figure5, plot1, plot2, plot3, plot4, residuals, Mouse_ID)
  
  