# Parasite community analysis
Field_par <- Field %>%
    filter(!is.na(predicted_weight_loss)) %>%
    mutate(
        infected_Aspiculuris = !is.na(Aspiculuris_sp) & Aspiculuris_sp > 0,
        infected_Syphacia = !is.na(Syphacia_sp) & Syphacia_sp > 0,
        infected_Trichuris = !is.na(Trichuris_muris) & Trichuris_muris > 0,
        infected_Mastophorus = !is.na(Mastophorus_muris) & Mastophorus_muris > 0
    ) %>%
    mutate(
        infected_Aspiculuris = as.factor(infected_Aspiculuris),
        infected_Syphacia = as.factor(infected_Syphacia),
        infected_Trichuris = as.factor(infected_Trichuris),
        infected_Mastophorus = as.factor(infected_Mastophorus)
    )

# Check sample sizes
cat("Sample size for parasite community analysis:", nrow(Field_par), "\n")
table(Field_par$infection_status)  # Use infection_status instead of MC.Eimeria

# Run the model
parasite_model <- lm(predicted_weight_loss ~ infection_status + infected_Aspiculuris 
                     + infected_Syphacia + infected_Trichuris + infected_Mastophorus, 
                     data = Field_par)

summary(parasite_model)
broom::tidy(parasite_model)

# Create coefficient plot
plot_summs(parasite_model, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "skyblue") 
ggsave(filename = paste0(an_fi, "/coefficient_plot_parasites.pdf"), 
       width = 8, height = 6, dpi = 300)
ggsave(filename = paste0(an_fi, "/coefficient_plot_parasites.jpeg"), 
       width = 8, height = 6, dpi = 300)
# Alternative: Create a summary from the model directly
parasite_model_summary <- summary(parasite_model)
parasite_clean_df <- as.data.frame(parasite_model_summary$coefficients)
parasite_clean_df$Term <- rownames(parasite_clean_df)
rownames(parasite_clean_df) <- NULL

# Clean up the term names
parasite_clean_df$Term <- case_when(
    parasite_clean_df$Term == "(Intercept)" ~ "Intercept",
    parasite_clean_df$Term == "infection_statusTRUE" ~ "Eimeria infection",
    parasite_clean_df$Term == "infected_AspiculurisFALSE" ~ "Aspiculuris sp.",
    parasite_clean_df$Term == "infected_SyphaciaTRUE" ~ "Syphacia sp.",
    parasite_clean_df$Term == "infected_TrichurisTRUE" ~ "Trichuris muris",
    parasite_clean_df$Term == "infected_MastophorusTRUE" ~ "Mastophorus muris"
)

# Reorder columns
parasite_final_table <- parasite_clean_df[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]

# Convert to gt table object first
parasite_gt_table <- parasite_supp_table %>%
    gt() %>%
    tab_header(
        title = "",
        subtitle = "Linear regression model results"
    ) %>%
    cols_label(
        Term = "Predictor",
        Estimate = "Estimate",
        `Std. Error` = "Std. Error", 
        `t-value` = "t-value",
        `p-value` = "p-value",
        `95% CI` = "95% CI"
    ) %>%
    fmt_number(
        columns = c("Estimate", "Std. Error", "t-value"),
        decimals = 3
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(rows = `p-value` == "< 0.001")
    )

# Now try the save function with the gt object
save_table_all_formats(parasite_gt_table, "Supplementary_Table_S11_parasite_community")

