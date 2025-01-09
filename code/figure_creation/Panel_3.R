importance_plot
predictions_random_for_lab
predictions
raincloud_plots__eimeria
int_MC


# Combine the plots

panel <-   
    (predictions_random_for_lab | importance_plot) /
    (coef_A | int_MC) /
    (Cor_infection_wl | Cor_OPG_wl ) +
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
       panel, width = 18, height =18, dpi = 300)




