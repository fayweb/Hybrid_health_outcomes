

## Let'S dive into it further
# what about a combination of the oocyst and infection intensity data with qpcr?
model <- lm(predicted_WL ~  infection_status * delta_ct_cewe_MminusE, data = Field)
summary(model)
confint(model)
modelsummary(model)

# Plot the summary with enhanced aesthetics
plot_summs(model, 
           scale = TRUE, 
           robust = TRUE,
           inner_ci_level = 0.95, 
           outer_ci_level = 0.99,
           coefs = c("Infected with Eimeria spp." = "infection_statusTRUE", 
                     "Infection intensity with Eimeria spp." = 
                         "infection_statusTRUE:delta_ct_cewe_MminusE")) +
    theme_minimal() +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
    ) +
    labs(
        x = "Estimate",
        y = "Predictor"
    ) +
    scale_color_manual(values = c("blue", "red")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") -> int_MC

int_MC

ggsave(paste0(an_fi, "/intensity_melting_curve.jpeg"), int_MC)


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
        x = "Eimeria spp. subspecies",
        y = "Predicted Weight Loss",
        color = "Species",
        shape = "Species"
    ) +
    theme_minimal() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_markdown(),
        legend.text = element_markdown()
    ) -> predictions


predictions



# Save the plot to a file
ggsave(paste0(an_fi, "/predicted_weight_loss_species.jpeg"),
       width = 8, height = 6, dpi = 300)


#### Are there any differences in prediced weight loss when we control 
# for other intestinal parasites?
# Eventhough our model is trained and tested on eimeria infections
Field <- Field %>%
    mutate(infection_intensity_Eim = delta_ct_cewe_MminusE)

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


##############remove eimeria
modelB <- lm(predicted_WL ~  infection_status + infected_Aspiculuris 
             + infected_syphasia + infected_crypto, data = Field_par)

summary(modelB)
modelsummary(modelB)

plot_summs(modelB, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "skyblue") -> coef_B
coef_B

ggsave(filename = paste0(an_fi, "/coefficient_plot_parasites.jpeg"), 
       width = 8, height = 6, dpi = 300)



#### just eimeria
modelc <- lm(predicted_WL ~ species_Eimeria * FEC_Eim_Ct, Field )
summary(modelc)

