#----------------------------------------------------------*
# Model 1: Infection intensity and predicted weight loss
model_intensity <- lm(predicted_WL ~ MC.Eimeria, data = Field)

# Model 2: Effect of infection species on predicted weight loss
model_species <- lm(predicted_WL ~ species_Eimeria, data = Field)

# Model 3: Effect of infection status and infection intensity
model_interaction <- lm(predicted_WL ~ infection_status * delta_ct_cewe_MminusE, data = Field)

# Create a tidy dataframe with model coefficients
model_list <- list("Infection Intensity" = model_intensity,
                   "Eimeria Species Effect" = model_species,
                   "Infection Status & Intensity Interaction" = model_interaction)

# Generate a summary table for Google Docs
modelsummary(model_list, 
             stars = TRUE, # Add significance stars
             output = "table_models.docx", # Save as Word file for Google Docs
             gof_map = c("r.squared", "adj.r.squared", "nobs"), # Include R² and sample size
             title = "Regression Model Summary for Predicted Weight Loss")

# Convert models to a tidy dataframe
model_df <- bind_rows(tidy(model_intensity) %>% mutate(Model = "Infection Intensity"),
                      tidy(model_species) %>% mutate(Model = "Eimeria Species Effect"),
                      tidy(model_interaction) %>% mutate(Model = "Infection Status & Intensity"))

# Save as CSV
write.csv(model_df, "output/tables/regression_models.csv", row.names = FALSE)
