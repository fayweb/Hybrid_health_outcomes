

## hybrid index + infectopm
Field$MC.Eimeria <- factor(Field$MC.Eimeria, levels = c("FALSE", "TRUE"))


################## Linearizing the HI
# calculating the expected heterozygocity
Field <- Field %>%
    mutate(HE = 2*HI*(1-HI), #linearize HI
           tolerance = predicted_WL / delta_ct_cewe_MminusE) 

model1 <- lm(predicted_WL ~ MC.Eimeria *delta_ct_cewe_MminusE * HE + 
                 HI + HE, Field)
summary(model1)

plot_summs(model1, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue") -> plot1

plot1

ggsave(filename = paste0(an_fi, "/coefficient_plot_model1.jpeg"), plot = plot1, 
       width = 5, height = 4, dpi = 300)



## hybrid index + infectopm
model2 <- lm(predicted_WL ~  HI + HE, Field)
summary(model2)
plot_summs(model2,  plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue")  -> plot_2

plot_2

ggsave(filename = paste0(an_fi, "/coefficient_plot_HE.jpeg"), plot = plot_2, 
       width = 5, height = 4, dpi = 300)

model3 <- lm(predicted_WL ~  MC.Eimeria, Field)
summary(model3)

plot_summs(model1, model2, model3, robust = TRUE, 
           scale = TRUE) -> model1_2
model1_2

ggsave(filename = paste0(an_fi, "figures/coefficient_plot_model1_2_3.jpeg"), 
       plot = model1_2, 
       width = 8, height = 6, dpi = 300)
#############################################
############################################
############################################

# Define colors
colors <- c("TRUE" = "forestgreen", "FALSE" = "purple")


Field %>%
    drop_na(MC.Eimeria)%>%
ggplot(aes(y = MC.Eimeria, x = predicted_WL, fill = MC.Eimeria)) + 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) + 
    scale_fill_manual(values = colors) +
    geom_boxplot(
        width = .15, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    stat_dots(
        # ploting on left side
        side = "left",
        # adjusting position
        justification = 1.1,
        # adjust grouping (binning) of observations
        binwidth = 0.25,
        alpha = 0.5) +
    geom_point(
        shape = 95,
        size = 15,
        alpha = .2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    coord_cartesian(ylim = c(1.2, 2.9), clip = "off") +
    theme_minimal() +
    labs(y = "Infection status with Eimerai spp.", 
         x = "Predicted weight loss" , 
         fill = "Infection status with Eimeria spp.") -> raincloud_plots__eimeria

raincloud_plots__eimeria

ggsave(plot = raincloud_plots__eimeria, filename = "figures/raincloud_eimeria.jpeg", 
       width = 6, 
       height = 4, dpi = 1000)


# Combine the plots
panel <- 
    (raincloud_plots__eimeria | model1_2) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel <- panel + 
    plot_annotation(title = 'Fig.', 
                    theme = theme(plot.title = element_text(size = 13, hjust = 0)))

# Display the panel figure
print(panel)

# Save the panel figure
ggsave('figure_panels/infected_hybrids.jpeg', 
       panel, width = 12, height = 5, dpi = 300)


