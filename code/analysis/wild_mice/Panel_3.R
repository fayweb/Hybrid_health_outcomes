# importance plot
# predictions correlations with observations in the lab
# predictions species eimeria
# 


importance_plot
predictions_random_for_lab
predictions
raincloud_plots__eimeria
int_MC


# Combine the plots

panel <-   
    (predictions_random_for_lab | importance_plot) /
    (predictions | raincloud_plots__eimeria | int_MC) +
    plot_annotation(tag_levels = 'A')


# Add a figure title
panel <- panel + 
    plot_annotation(title = 'Fig. 6', 
                    theme = theme(plot.title = 
                                      element_text(size = 13, hjust = 0)))

# Display the panel figure
print(panel)

# Save the panel figure
ggsave(paste0(panels_fi, '/infection_predictions_panels.jpeg'), 
       panel, width = 16, height =10, dpi = 300)




