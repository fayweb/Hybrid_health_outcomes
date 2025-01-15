# 7.7: Testing impact of infection status with Eimeria spp. on wild mice
# Do we predict higher weight loss for mice that are infected with Eimeria spp? 
# Is there any evidence for an association between infection status, hybridicity
# and their impact on weight loss? 

# make infection status a factor
#Field$MC.Eimeria <- factor(Field$MC.Eimeria, levels = c("FALSE", "TRUE"))
#############################################
############################################
############################################
## How does the predicted weight loss differ betwen infected species
species <- lm(predicted_WL ~ species_Eimeria, data = Field)
summary(species)
modelsummary(species)


# Create predicted values using ggpredict
preds <- ggpredict(species, terms = "species_Eimeria")


# Plot the predicted values
ggplot(preds, aes(x = x, y = predicted, color = x)) +
    geom_point(#aes(shape = x), 
        size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, size = 0.7) +
    geom_line(aes(group = group, color = "black")) +
    scale_color_manual(values = color_mapping_f, labels = labels_f) +
    labs(
        # title = "Predicted weight loss in wild mice, 
        #categorized by the species of Eimeria",
        x = "Eimeria spp. subspecies",
        y = "Predicted Weight Loss",
        color = "Species",
        shape = "Species"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_markdown(),
        legend.text = element_markdown()
    ) -> predictions


predictions



# Save the plot to a file
ggsave(paste0(an_fi, "/predicted_weight_loss_species.jpeg"),
       #width = 8, height = 6, 
       dpi = 1000)


#### Are there any differences in prediced weight loss when we control 
# for other intestinal parasites?
# Eventhough our model is trained and tested on eimeria infections

Field_par <- Field %>%
    mutate(
        infected_Aspiculuris = Aspiculuris_sp != 0,
        infected_syphasia = Syphacia_sp != 0,
        infected_crypto = ILWE_Crypto_Ct != 0
    )

Field_par <- Field_par %>%
    mutate(
        infected_Aspiculuris = as.factor(infected_Aspiculuris),
        infected_syphasia = as.factor(infected_syphasia),
        infected_crypto = as.factor(infected_crypto)
    )

modelA <- lm(predicted_WL ~ MC.Eimeria + infected_Aspiculuris 
             + infected_syphasia + infected_crypto, data = Field_par)

summary(modelA)
modelsummary(modelA)

plot_summs(modelA, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "skyblue") -> coef_A
coef_A

ggsave(filename = paste0(an_fi, "/coefficient_plot_parasites.jpeg"), 
       width = 8, height = 6, dpi = 300)


###### Plotting the raincloud plot to show the effect of Eimeria spp.
# on predicted weight loss
# Define colors
gsub(pattern = "TRUE",replacement =  "Infected with Eimeria spp.", 
     x = Field$MC.Eimeria) -> Field$MC.Eimeria
gsub(pattern = "FALSE",replacement =  "Uninfected", 
     x = Field$MC.Eimeria) -> Field$MC.Eimeria



colors <- c("Infected with Eimeria spp." = "red", "Uninfected" = "steelblue")


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
    labs(
        #title = "Mice infected with Eimeria spp. experience 
        #significantly greater predicted weight loss",
        y = "Infection status with Eimerai spp.", 
        x = "Predicted weight loss" , 
        fill = "Infection status with Eimeria spp.") +
    theme(legend.position = "blank") -> raincloud_plots__eimeria

raincloud_plots__eimeria

ggsave(plot = raincloud_plots__eimeria, filename = 
           paste0(an_fi, "/raincloud_eimeria.jpeg"), 
       width = 6, 
       height = 4, 
       dpi = 1000)


