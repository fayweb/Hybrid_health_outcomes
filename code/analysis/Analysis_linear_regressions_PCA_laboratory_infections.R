# 6.2: PCA 
# Regressions with pc axes 
# Plots: pc1_current_infection, pc2_current_infection, coefs5
#----------------------------------------------------------*
model_1 <- lm(WL_max ~ PC1 + PC2 +
                  current_infection + delta_ct_cewe_MminusE +
                  mouse_strain + immunization + 
                  weight_dpi0, data = lab )

summary(model_1)

tab_model(model_1)

stargazer(model_1,
          type = "text", out = paste0(tables, 
                                      "/predictors_weightloss.doc"), 
          title = "Linear models - Predicting maximum weight loss")


# Extract the residuals from the model
residuals <- resid(model_1)

# Create a data frame with the residuals
residuals_df <- data.frame(residuals = residuals)

# Create the QQ plot
residuals_1 <-
    ggplot(residuals_df, aes(sample = residuals)) +
    stat_qq(color = "blue") +
    ggtitle("QQ Plot of Residuals") +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    theme_minimal()

residuals_1

ggsave(filename = paste0(an_fi, "/residuals_model_1.jpeg"), 
       plot = residuals_1, 
       width = 12, height = 6, dpi = 600)

# Extract the fitted values from the model
fitted_values <- fitted(model_1)

# Create a data frame with the residuals and the fitted valueshttp://127.0.0.1:11745/graphics/plot_zoom_png?width=674&height=334
data_df <- data.frame(residuals = residuals, fitted_values = fitted_values)

# Create the scatter plot
residuals_vs_fitted <-
    ggplot(data_df, aes(x = fitted_values, y = residuals)) +
    geom_point(color = "blue") +
    ggtitle("Residuals vs Fitted Values") +
    xlab("Fitted Values") +
    ylab("Residuals") +
    theme_minimal()

residuals_vs_fitted

ggsave(filename = paste0(an_fi, "/residuals_vs_fitted.jpeg"), 
       plot = residuals_vs_fitted, 
       width = 12, height = 6, dpi = 600)


# without host data
model_2 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE, data = lab)

summary(model_2)

# only pc1 + pc2
model_3 <- lm(WL_max ~ PC1 + PC2 , data = lab)

summary(model_3)


model_4 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE + 
                  weight_dpi0, data = lab)


summary(model_4)
plot_coefs(model_1, model_2, model_3) -> coef_plot_pca
coef_plot_pca


## make a coefficient plot with just the PC1 and PC2
coef_plot_PC1PC2 <- plot_coefs(model_3, plot.distributions = TRUE,
                               rescale.distributions = TRUE,
                               colors = "blue")
coef_plot_PC1PC2

ggsave(filename = paste0(an_fi, "/plot_sums_pc1_pc2_simple.pdf"),
       plot = coef_plot_PC1PC2, width = 6, height = 4, dpi = 300)

# remove gene information
model_5 <- lm(WL_max ~  current_infection + delta_ct_cewe_MminusE +
                  mouse_strain + immunization + 
                  weight_dpi0, data = lab )

summary(model_5)

## Please cite as:
##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
stargazer(model_1, model_5, model_3,
          type = "text",
          out = paste0(tables, "/stargazer.docx"), 
          title = "Linear models - Predicting maximum weight loss",
          align = TRUE)

export_summs(model_1, model_5, model_3,
             scale = TRUE, to.file = "docx", 
             file.name = paste0(tables, "/lab_model1_3.docx"))

models <- list(
    "Model 1" = model_1,
    #"Model 2" = model_2,
    "Model 2" = model_5,
    "Model 3" = model_3)

summ(model_1)
modelsummary(models = models) -> s_models
str(s_models)


modelsummary(models, stars = TRUE, #gof_omit = "IC|Adj|F|RMSE|Log", 
             output = paste0(tables, "/lab_model1_3.docx"))

modelsummary(models = as.list(model_1, model_5, model_3), 
             output = paste0(tables, "/lab_model1.docx"))

summ(model_1)
#modelsummary(models = c(model_1, model_5, model_3))

plot_summs(model_1, model_5, model_3, 
           colors = "CUD", rescale.distributions = TRUE, 
           plot.distributions = FALSE, point.size = 3, 
           robust = TRUE,
           model.names = c("Complex model", "Without PC1 and PC2", 
                           "Solely PC1 and PC2")) -> coefs1_3

coefs1_3

ggsave(filename = paste0(an_fi, "/plot_sums_complex_pca.jpeg"),
       plot = coefs1_3, width = 6, height = 4)

####################################
# Tidy the models
tidy_model_1 <- tidy(model_1)
tidy_model_5 <- tidy(model_5)
tidy_model_3 <- tidy(model_3)

# Combine the tidied models with model labels
all_models <- bind_rows(
    tidy_model_1 %>% mutate(model = "Model 1"),
    tidy_model_5 %>% mutate(model = "Model 2"),
    tidy_model_3 %>% mutate(model = "Model 3")
)

# Group mouse_strain terms into a single summary row per model
mouse_strain_summary <- all_models %>%
    filter(grepl("mouse_strain", term)) %>%
    group_by(model) %>%
    summarise(
        term = "Mouse Strain (combined)",
        estimate = paste0(round(estimate, 3), collapse = ", "),
        std.error = paste0("(", round(std.error, 3), ")", collapse = ", "),
        p.value = paste0(round(p.value, 3), collapse = ", ")
    ) %>%
    mutate(across(c(estimate, std.error, p.value), as.character))  # Ensure consistency in data type

# Filter out the mouse_strain terms from the non-strain rows
non_mouse_strain_terms <- all_models %>%
    filter(!grepl("mouse_strain", term)) %>%
    mutate(across(c(estimate, std.error, p.value), as.character))  # Convert to character to match

# Combine non-mouse_strain terms with mouse_strain_summary
final_table <- bind_rows(non_mouse_strain_terms, mouse_strain_summary)

# Arrange rows by model and move Mouse Strain row down
final_table <- final_table %>%
    arrange(model, desc(term != "Mouse Strain (combined)"))

# Create a final clean table
publication_table <- final_table %>%
    dplyr::select(model, term, estimate, std.error, p.value) %>%  # Select relevant columns
    kbl(format = "html", digits = 3, caption = "Compact Regression Table with Mouse Strain Combined") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
    add_header_above(c(" " = 2, "Model 1" = 1, "Model 2" = 1, "Model 3" = 1))  # Adjusted for 5 total columns

# Save as HTML
save_kable(publication_table, file = "output/tables/compact_regression_table_combined.html")

# Convert to PNG for Google Docs
webshot("output/tables/compact_regression_table_combined.html", 
        file = "output/tables/compact_regression_table_combined.png")


#########################################################################
##########################################################################
########################################################################

model_6 <- lm(WL_max ~ PC1 * current_infection + PC2 *current_infection, 
              data = lab)

summary(model_6)

## Please cite as:
##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.


# Create a regression table using stargazer
stargazer(model_6, type = "text",
          title = "Table X: Regression Results for Weight Loss Prediction",
          dep.var.labels = "Maximum Weight Loss (WL_max)",
          covariate.labels = c("PC1", "E. ferrisi", "E. falciformis", "PC2", 
                               "PC1 × E. ferrisi", "PC1 × E. falciformis", 
                               "PC2 × E. ferrisi", "PC2 × E. falciformis"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "Significance levels: *** p<0.001, ** p<0.01, * p<0.05"
)


stargazer(model_6, type = "latex", out = "table_x.tex",
          title = "Regression Results for Weight Loss Prediction",
          dep.var.labels = "Maximum Weight Loss (WL_max)",
          covariate.labels = c("PC1", "E. ferrisi", "E. falciformis", "PC2", 
                               "PC1 × E. ferrisi", "PC1 × E. falciformis", 
                               "PC2 × E. ferrisi", "PC2 × E. falciformis"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "Significance levels: *** p<0.001, ** p<0.01, * p<0.05"
)



# Extract coefficients, standard errors, and p-values
results_df <- tidy(model_6)

# Add significance stars
results_df$Significance <- ifelse(results_df$p.value < 0.001, "***",
                                  ifelse(results_df$p.value < 0.01, "**",
                                         ifelse(results_df$p.value < 0.05, "*", "")))

# Rename columns for clarity
colnames(results_df) <- c("Variable", "Estimate", "Std. Error", "t value", "p-value", "Significance")

# Save as CSV for Google Docs
write.csv(results_df, "output/tables/table_x.csv", row.names = FALSE)


plot_summs(model_6) -> coefs6

coefs6

ggsave(filename = paste0(an_fi, "/plot_sums_mix_PCA.jpeg"),
       plot = coefs6, width = 6, height = 4)
#see the ggefects
effects <- ggpredict(model_6)

pc1_current_infection <- 
    ggpredict(model_6, terms = c("PC1")) %>% 
    plot(colors = "darkorchid") +   # Use a refined shade of blue
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    #  ggtitle("Effect of PC1 on Predicted Weight Loss") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )


pc1_current_infection

ggsave(filename = paste0(an_fi, "/pc1_current_infection.jpeg"), 
       plot = pc1_current_infection, 
       width = 6, height = 4, dpi = 1000)

pc2_current_infection <- 
    ggpredict(model_6, terms = c("PC2")) %>% 
    plot(colors = "darkorchid") +   # Use a refined shade of blue
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )

pc2_current_infection

ggsave(filename = paste0(an_fi, "/pc2_current_infection.jpeg"), 
       plot = pc2_current_infection, 
       width = 6, height = 4, dpi = 1000)


plot_summs(model_6) -> coef_interaction
modelsummary(model_6, stars = TRUE, gof_omit = "IC|Adj|F|RMSE|Log", 
             output = paste0(tables, "/mixed_effects_pca.docx"))

# Now create the scatter plot using this color mapping
ggpredict(model_6, terms = c("PC1", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    labs(color = "Infection group") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()) -> pc1_WL_current_infection

pc1_WL_current_infection

ggsave(paste0(an_fi, "/pc1_WL_current_infection.jpeg"), pc1_WL_current_infection, 
       width = 8, height = 6, dpi = 1000)


# Now create the scatter plot using this color mapping
# Now create the scatter plot using this color mapping
ggpredict(model_6, terms = c("PC2", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(color = "Infection group") +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()) -> pc2_WL_current_infection

pc2_WL_current_infection

ggsave(paste0(tables, "/pc2_WL_current_infection.jpeg"),
       pc2_WL_current_infection, width = 8, height = 6, dpi = 1000)

####################
################### Create the panel figure
figure_panel <- ggarrange(pca_variables, biplot, 
                          pc1_current_infection, pc2_current_infection,
                          pc1_WL_current_infection, pc2_WL_current_infection,
                          labels = c("A", "B", "C", "D", "E", "F"),
                          ncol = 2, nrow = 3)

# Adding the title "Figure 1" to the entire arrangement
figure_panel <- annotate_figure(figure_panel, 
                                top = text_grob("Fig. 2", size = 14, 
                                                face = "bold"))


ggsave(paste0(panels_fi, "/panel_regression_pca.jpeg"), 
       figure_panel, width = 12, height = 10, dpi = 300)


###########################biplot and linear models
panel_biplot_regr <- ggarrange(biplot, coefs1_3,
                               labels = c("A", "B"),
                               ncol = 2)

panel_biplot_regr <- annotate_figure(panel_biplot_regr,
                                     top = text_grob("Fig. 3", size = 14, 
                                               face = "bold"))

ggsave(paste0(panels_fi, "/panel_regression_biplot.jpeg"), 
       panel_biplot_regr, width = 21, height = 10, dpi = 300)

################### Create the simplified figure# combine
panel_figure5 <- 
    #(pc1_current_infection | pc2_current_infection ) /
    (pc1_WL_current_infection | pc2_WL_current_infection) /
    free(coefs6) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure5 <- panel_figure5 + 
    plot_annotation(title = 'Fig. 4', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))# +
    #plot_layout(heights = c(1, 1,1), 
     #           widths = c(1,1,1))

# Save the panel figure
ggsave(paste0(panels_fi, "/panel_regression_pca_interaction.jpeg"), 
       panel_figure5, width = 11, height = 8, dpi = 300)

### Create the panel figure
PCA_panel <-
    (biplot | coefs1_3 ) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
PCA_panel <- PCA_panel + 
    plot_annotation(title = 'Fig. 3', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))



# Save the panel figure
ggsave(paste0(panels_fi, "/panel_pca.jpeg"), 
       PCA_panel, width = 18, height = 8, dpi = 300)

#rm(residuals_1, residuals_df, residuals_vs_fitted, model_1, model_2,
 #  model_3, model_4, model_5, models, effects, data_df)
#rm(circ, mouse, pca.vars, pca.vars.m, 
 #  pca_var, var.contrib.matrix, res.pca, 
  # var.contrib, pca_variables, coef6, coef_interaction, contributions_pc1,
   #contributions_pc2)

