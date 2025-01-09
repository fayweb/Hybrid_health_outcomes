# Laboratory Coefficient plot
lab_plot <- coef_mmr
lab_plot


# Field coefficient plot
field_plot <- coef_mmr_B
field_plot


pca_plot


#######################
# combine
comb <- ((lab_plot | field_plot) /
             pca_plot) +
    #  plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
comb <- comb + 
    plot_annotation(title = 'Figure 2', 
                    theme = theme(plot.title = element_text(size = 16, hjust = 0)))



# Display the panel figure
# print(comb)

# Save the panel figure
ggsave(paste0(panels_fi, "/panel_bridging_immune_expression.jpeg"), 
       comb, width = 12, height = 10, dpi = 300)    
