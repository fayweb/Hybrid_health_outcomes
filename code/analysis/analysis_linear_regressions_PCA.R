###PC1 PC2 linear regression
lab$current_infection <- factor(lab$current_infection, 
                                levels = c("uninfected", "E_falciformis", "E_ferrisi"))

lab$mouse_strain <- as.factor(lab$mouse_strain)
lab$immunization <- factor(lab$immunization, levels = 
                               c("naive", "uninfected",
                                 "heterologous", "homologous"))

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
plot_coefs(model_1, model_2, model_3, model_4)

# remove gene information
model_5 <- lm(WL_max ~  current_infection + delta_ct_cewe_MminusE +
                  mouse_strain + immunization + 
                  weight_dpi0, data = lab )

summary(model_5)

## Please cite as:
##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
stargazer(model_1, model_2, model_3, model_5,
          type = "text",
          out = paste0(tables, "/stargazer.docx"), 
          title = "Linear models - Predicting maximum weight loss",
          align = TRUE)

export_summs(model_1, model_2, model_3,model_5,
             scale = TRUE, to.file = "docx", 
             file.name = paste0(tables, "/lab_model1_3.docx"))

models <- list(
    "Model 1" = model_1,
    "Model 2" = model_2,
    "Model 3" = model_3,
    "Model 4" = model_5)

modelsummary(models, stars = TRUE, gof_omit = "IC|Adj|F|RMSE|Log", 
             output = paste0(tables, "/lab_model1_3.docx"))

modelsummary(models = as.list(model_1, model_2, model_3), 
             output = paste0(tables, "/lab_model1.docx"))


model_6 <- lm(WL_max ~ PC1 * current_infection + PC2 *current_infection, 
              data = lab)

summary(model_6)
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
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) -> pc1_WL_current_infection

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
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) -> pc2_WL_current_infection

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

################### Create the simplified figure# combine
panel_figure5 <- 
    (pc1_current_infection | pc2_current_infection ) /
    (pc1_WL_current_infection | pc2_WL_current_infection) /
    free(coef_interaction) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure5 <- panel_figure5 + 
    plot_annotation(title = 'Fig. 5', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0))) +
    plot_layout(heights = c(1, 1,1), 
                widths = c(1,1,1))

# Save the panel figure
ggsave(paste0(panels_fi, "/panel_regression_pca_interaction.jpeg"), 
       panel_figure5, width = 13, height = 12, dpi = 300)



rm(residuals_1, residuals_df, residuals_vs_fitted, model_1, model_2,
  model_3, model_4, model_5, models, effects, data_df)
rm(circ, mouse, pca.vars, pca.vars.m, 
   pca_var, var.contrib.matrix, res.pca, 
   var.contrib, pca_variables, coef6, coef_interaction, contributions_pc1,
   contributions_pc2)
