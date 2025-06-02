# ***********************************************************
# Script 3: random_forest_lab_validation.R
# ***********************************************************

# ***********************************************************
# Random Forest Laboratory Validation
# Purpose: Test RF predictions against infection parameters in lab data
# Author: Fay Webster
# ***********************************************************

cat("=== RANDOM FOREST LABORATORY VALIDATION ===\n")
## ----predicting_weight_loss_model---
set.seed(333)

#train the model
WL_predict_gene <- randomForest(WL_max ~., data = train.data, 
                                proximity = TRUE, ntree = 26) 

# ntree = number of trees     
# save the model 
saveRDS(WL_predict_gene, file =  paste0(cmodels, "WL_predict_gene_training_data_set.RDS"))
print(WL_predict_gene)

predict_WL_cv <- rf.crossValidation(x = WL_predict_gene, xdata = train.data, 
                                    p = 0.10, n = 99, ntree = 26)
#The predict() function in R is used to predict the values based on the 
# input data.
predictions <- predict(WL_predict_gene, test.data)

# assign test.data to a new object, so that we can make changes
result <- test.data

#add the new variable of predictions to the result object
result <- cbind(result, predictions)

# what is the correlation between predicted and actual data?

cor.test(result$WL_max, result$predictions)


test_lab <- lab %>%
    left_join(result, by = c("WL_max", "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                             "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                             "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                             "TICAM1", "TNF"))

test_lab <- test_lab %>%
    drop_na(predictions)



#######################################
##############################################
#####################################

model <- lm(predictions ~ WL_max * current_infection, data = test_lab)    
summary(model)

#### Plotting
ggpredict(model, terms = c("WL_max", "current_infection")) %>% 
    plot(colors = "darkorchid") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss during infections") +
    theme_minimal() +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.2),
        legend.text = element_markdown()) -> lm_weight_loss_predictions

lm_weight_loss_predictions

ggsave(filename = paste0(an_fi, "/obs_pred_treatment_groups_random.pdf"),
       width = 6, height = 5, dpi = 300)

############################
##############uusing only the current infection to predict predicted weight 
# loss
model_c <- lm(predictions ~ current_infection, data = test_lab)    
summary(model_c)
preds <- ggpredict(model_c, terms = "current_infection")

#### Plotting
ggplot(preds, aes(x = x, y = predicted, color = x)) +
    geom_point(#aes(shape = x), 
        size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, size = 0.7) +
    geom_line(aes(group = group, color = "black")) +
    scale_color_manual(values = color_mapping, labels = labels) +
    labs(
        # title = "Bar plot showing the predicted maximum weight 
        #loss for each infection group",
        x = "Experimental infection groups",
        y = "Predicted maximum weight loss",
        color = "current_infection",
        shape = "current_infection"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_markdown(),
        legend.text = element_markdown()
    ) -> lm_weight_loss_predictions_c

lm_weight_loss_predictions_c

ggsave(filename = paste0(an_fi, "/random_foreest_lab_predictions_eimeria.pdf"),
       width = 6, height = 5, dpi = 300)


#######################################
#####################################
########################################
model <- lm(predictions ~ WL_max * delta_ct_cewe_MminusE , data = test_lab %>%
                filter(MC.Eimeria == "TRUE"))
summary(model)

ggpredict(model, terms = c("WL_max", "delta_ct_cewe_MminusE"), interactive=TRUE) %>% 
    plot() +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss") +
    theme_minimal() +
    scale_color_manual(values = c("darkred", "gold", "violet")) +
    scale_fill_manual(values = c("darkred", "gold", "violet")) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))-> lm_short

lm_short

ggsave(filename = paste0(an_fi, "/deltact_random.pdf"),
       width = 6, height = 5, dpi = 300)


#######################################
#####################################
########################################
model_d <- lm(predictions ~  delta_ct_cewe_MminusE , data = test_lab %>%
                  filter(MC.Eimeria == "TRUE"))
summary(model_d)

ggpredict(model_d, terms = c("delta_ct_cewe_MminusE"), interactive=TRUE) %>% 
    plot(color = "purple") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Relationship between infection intensity 
    #         and predicted maximum weight loss") +
    xlab("Infection intensity with Eimeria delta Ct") +
    ylab("Predicted maximum weight loss") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))-> lm_short

lm_short

preds <- ggpredict(model_d, terms = "delta_ct_cewe_MminusE")
ggsave(filename = paste0(an_fi, "/deltact_random.pdf"),
       width = 6, height = 5, dpi = 300)
#### Plotting
ggplot(preds, aes(x = x, y = predicted, color = x)) +
    # geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, size = 0.7) +
    geom_abline() +
    # scale_color_manual(values = color_mapping, labels = labels) +
    labs(
        x = "Infection group",
        y = "Predicted Weight Loss",
        color = "current_infection",
        shape = "current_infection"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_markdown(),
        legend.text = element_markdown()
    )

ggsave(filename = paste0(an_fi, "/deltact_random.pdf"),
       width = 6, height = 5, dpi = 300)


summary(model_d)
### plotting
test_lab %>%
    ggplot(aes(x = predictions, y = WL_max, color = current_infection)) +
    geom_point(aes(size = 5),#delta_ct_cewe_MminusE), 
               alpha = 0.7) +
    labs(
        x = "Predictions: Maximum weight loss", 
        y = "Observed: Maximum weight loss",
        #  title = "Relationship between Predicted and Observed Weight Loss",
        #subtitle = "Grouped by Current Infection and Sized by Delta CT Value",
        color = "Parasite strain",
        # size = "Caecal infection intensities, Delta Ct value",
        #shape = "Delta Ct treshold"
    ) +
    theme_minimal() +
    theme(
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_size_continuous(range = c(2, 10)) +
    guides(size = "none") +
    theme(legend.position = c(0.8, 0.2),
          legend.text = element_markdown())-> predictions_random_for_lab

predictions_random_for_lab


ggsave(plot = predictions_random_for_lab, 
       filename = paste0(an_fi, "/predictions_random_for_lab.pdf"), 
       width = 6, height = 5,
       dpi = 1000)


combi_plot <- (importance_plot | predictions_random_for_lab) /
    (lm_weight_loss_predictions_c | lm_short) +
    #plot_layout(guides = 'collect') + # Collect all legends into a 
    #single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

combi_plot

# Add a figure title
combi_plot <- combi_plot + 
    plot_annotation(title = 'Fig. 5', 
                    theme = theme(plot.title = element_text(size = 13, 
                                                            hjust = 0)))

# Display the panel figure
print(combi_plot)


ggsave(plot = combi_plot, 
       filename = paste0(panels_fi, "/variableimp_rand_results_lab.pdf"), 
       width = 12, 
       height = 8, dpi = 1000)

# Calculate the linear model
lm_fit <- lm(WL_max ~ predictions, data = test_lab)

# Extract coefficients for the model formula
intercept <- round(coef(lm_fit)[1], 2)
slope <- round(coef(lm_fit)[2], 2)
formula_text <- paste0("WL_max = ", intercept, " ", 
                       ifelse(slope >= 0, "+ ", "- "), 
                       abs(slope), " * predictions")

# Calculate correlation
cor_value <- round(cor(test_lab$WL_max, test_lab$predictions), 2)
cor_text <- paste0("Rho = ", cor_value)

test_lab   %>%
    ggplot(aes(x = predictions, y = WL_max)) +
    geom_smooth(method = lm, se = TRUE) +
    labs(x = "Predictions: Maximum weight loss", 
         y = "Observed: Maximum weight loss") +
    geom_point(aes(x = predictions, y = WL_max, size = 0.8, alpha = 0.3)) +
    labs(x = "Predictions: Maximum weight loss", 
         y = "Observed: Maximum weight loss") +
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "none") +
    annotate("text", 
             x = min(test_lab$predictions), y = max(test_lab$WL_max), 
             label = formula_text, hjust = 0, 
             vjust = 4, 
             size = 4, color = "blue") +
    annotate("text", x = min(test_lab$predictions), 
             y = max(test_lab$WL_max), 
             label = cor_text, hjust = 0, vjust = 1.5, 
             size = 4, color = "blue") -> linear_plot

linear_plot

ggsave(filename = paste0(an_fi, "/linear_model_of_random_forest.pdf"), plot = linear_plot, 
       width = 10, height = 6,
       dpi = 1000)




figure_panel_2 <- ggarrange(predictions_random_for_lab,
                            ggarrange(importance_plot, lm_short, 
                                      labels = c("B", "C"), ncol = 2),
                            nrow = 2, labels = "A")



# Adding the title "Figure 1" to the entire arrangement
figure_panel_2 <- annotate_figure(figure_panel_2, 
                                  top = text_grob("Figure 2", size = 14, 
                                                  face = "bold"))

print(figure_panel_2)


ggsave(paste0(panels_fi, "/panel_random_forest_lab_alternative.pdf"), 
       figure_panel_2, 
       width = 18, height = 18, dpi = 300)

# ***********************************************************
# Create Supplementary Figure S2: Laboratory Validation
# ***********************************************************

# Panel A: Predicted weight loss by infection status (cleaner version)
validation_panel_A <- ggplot(test_lab, aes(x = current_infection, y = predictions)) +
    geom_boxplot(alpha = 0.7, aes(fill = current_infection)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(x = "Infection Status", 
         y = "Predicted Weight Loss (%)",
         title = "A) Validation: Infection Status") +
    theme_minimal() +
    theme(axis.text.x = element_markdown(),
          legend.position = "none")

# Panel B: Infection intensity correlation (cleaner version)  
infected_test <- test_lab %>% filter(MC.Eimeria == "TRUE", !is.na(delta_ct_cewe_MminusE))

validation_panel_B <- ggplot(infected_test, aes(x = delta_ct_cewe_MminusE, y = predictions)) +
    geom_point(alpha = 0.7, size = 2, color = "purple") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(x = "Infection Intensity (ΔCt)", 
         y = "Predicted Weight Loss (%)",
         title = "B) Validation: Infection Intensity") +
    annotate("text", 
             x = min(infected_test$delta_ct_cewe_MminusE) + 1,
             y = max(infected_test$predictions) * 0.9,
             label = paste("r =", round(cor(infected_test$delta_ct_cewe_MminusE, 
                                            infected_test$predictions, use = "complete.obs"), 3)),
             size = 4) +
    theme_minimal()

# Panel C: Test set performance summary
model_summary_data <- data.frame(
    Metric = c("Test Correlation", "Species Effect (E. ferrisi)", "Species Effect (E. falciformis)", "Intensity Correlation"),
    Value = c(cor(test_lab$WL_max, test_lab$predictions),
              coef(model_c)["current_infectionE. ferrisi"],
              coef(model_c)["current_infectionE. falciformis"], 
              cor(infected_test$delta_ct_cewe_MminusE, infected_test$predictions, use = "complete.obs")),
    P_value = c(cor.test(test_lab$WL_max, test_lab$predictions)$p.value,
                summary(model_c)$coefficients["current_infectionE. ferrisi", "Pr(>|t|)"],
                summary(model_c)$coefficients["current_infectionE. falciformis", "Pr(>|t|)"],
                cor.test(infected_test$delta_ct_cewe_MminusE, infected_test$predictions)$p.value)
)

validation_panel_C <- ggplot(model_summary_data, aes(x = Metric, y = Value)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    labs(x = "Validation Metric", 
         y = "Effect Size / Correlation",
         title = "C) Validation Summary") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine into Supplementary Figure S2
supp_fig_S2 <- grid.arrange(
    validation_panel_A, validation_panel_B, validation_panel_C,
    ncol = 2, nrow = 2,
    top = "Supplementary Figure S2: Random Forest Laboratory Validation"
)

# ***********************************************************
# Create Supplementary Figure S2: Laboratory Validation
# ***********************************************************

# Panel A: Predicted weight loss by infection status (cleaner version)
validation_panel_A <- ggplot(test_lab, aes(x = current_infection, y = predictions)) +
    geom_boxplot(alpha = 0.7, aes(fill = current_infection)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(x = "Infection Status", 
         y = "Predicted Weight Loss (%)",
         title = "A) Validation: Infection Status") +
    theme_minimal() +
    theme(axis.text.x = element_markdown(),
          legend.position = "none")

# Panel B: Infection intensity correlation (cleaner version)  
infected_test <- test_lab %>% filter(MC.Eimeria == "TRUE", !is.na(delta_ct_cewe_MminusE))

validation_panel_B <- ggplot(infected_test, aes(x = delta_ct_cewe_MminusE, y = predictions)) +
    geom_point(alpha = 0.7, size = 2, color = "purple") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(x = "Infection Intensity (ΔCt)", 
         y = "Predicted Weight Loss (%)",
         title = "B) Validation: Infection Intensity") +
    annotate("text", 
             x = min(infected_test$delta_ct_cewe_MminusE) + 1,
             y = max(infected_test$predictions) * 0.9,
             label = paste("r =", round(cor(infected_test$delta_ct_cewe_MminusE, 
                                            infected_test$predictions, use = "complete.obs"), 3)),
             size = 4) +
    theme_minimal()

# Panel C: Test set performance summary
model_summary_data <- data.frame(
    Metric = c("Test Correlation", "Species Effect (E. ferrisi)", "Species Effect (E. falciformis)", "Intensity Correlation"),
    Value = c(cor(test_lab$WL_max, test_lab$predictions),
              coef(model_c)["current_infectionE. ferrisi"],
              coef(model_c)["current_infectionE. falciformis"], 
              cor(infected_test$delta_ct_cewe_MminusE, infected_test$predictions, use = "complete.obs")),
    P_value = c(cor.test(test_lab$WL_max, test_lab$predictions)$p.value,
                summary(model_c)$coefficients["current_infectionE. ferrisi", "Pr(>|t|)"],
                summary(model_c)$coefficients["current_infectionE. falciformis", "Pr(>|t|)"],
                cor.test(infected_test$delta_ct_cewe_MminusE, infected_test$predictions)$p.value)
)

validation_panel_C <- ggplot(model_summary_data, aes(x = Metric, y = Value)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    labs(x = "Validation Metric", 
         y = "Effect Size / Correlation",
         title = "C) Validation Summary") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine into Supplementary Figure S2
supp_fig_S2 <- grid.arrange(
    validation_panel_A, validation_panel_B, validation_panel_C,
    ncol = 2, nrow = 2,
    top = "Supplementary Figure S2: Random Forest Laboratory Validation"
)


# Save Supplementary Figure S2
ggsave(filename = paste0(an_fi, "/supplementary_figure_S2_lab_validation.pdf"),
       plot = supp_fig_S2, width = 12, height = 8, dpi = 300)

# Also save individual panels
ggsave(filename = paste0(an_fi, "/supp_S2A_infection_status.pdf"),
       plot = validation_panel_A, width = 6, height = 5, dpi = 300)

ggsave(filename = paste0(an_fi, "/supp_S2B_infection_intensity.pdf"),
       plot = validation_panel_B, width = 6, height = 5, dpi = 300)
                         
# ***********************************************************
# Create Supplementary Table S2 using your function
# ***********************************************************
# Create validation statistics table
validation_stats_data <- data.frame(
    Analysis = c("Test Set Correlation", "E. ferrisi Effect", "E. falciformis Effect", 
                 "Infection Intensity", "Overall Model R²"),
    Estimate = round(c(cor(test_lab$WL_max, test_lab$predictions),
                       coef(model_c)["current_infectionE. ferrisi"],
                       coef(model_c)["current_infectionE. falciformis"],
                       coef(model_d)["delta_ct_cewe_MminusE"],
                       summary(model_c)$r.squared), 3),
    P_Value = round(c(cor.test(test_lab$WL_max, test_lab$predictions)$p.value,
                      summary(model_c)$coefficients["current_infectionE. ferrisi", "Pr(>|t|)"],
                      summary(model_c)$coefficients["current_infectionE. falciformis", "Pr(>|t|)"],
                      summary(model_d)$coefficients["delta_ct_cewe_MminusE", "Pr(>|t|)"],
                      NA), 4),
    n = c(nrow(test_lab), nrow(test_lab), nrow(test_lab), nrow(infected_test), nrow(test_lab))
)

# Create beautiful gt table
supp_table_S2 <- validation_stats_data %>%
    gt() %>%
    tab_header(
        title = "",
        subtitle = "Random Forest laboratory validation statistics"
    ) %>%
    tab_spanner(
        label = "Statistical Results",
        columns = c(Estimate, P_Value)
    ) %>%
    cols_label(
        Analysis = "Validation Analysis",
        Estimate = "Estimate/Correlation", 
        P_Value = "P-value",
        n = "Sample Size"
    ) %>%
    fmt_number(
        columns = c(Estimate),
        decimals = 3
    ) %>%
    fmt_scientific(
        columns = c(P_Value),
        decimals = 2
    ) %>%
    tab_source_note(
        source_note = "Test set validation results for random forest model predictions against infection parameters"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>%
    tab_style(
        style = cell_fill(color = "lightgray"),
        locations = cells_body(rows = P_Value < 0.05, columns = P_Value)
    ) %>%
    # Make Eimeria species names italic
    tab_style(
        style = cell_text(style = "italic"),
        locations = cells_body(
            columns = Analysis,
            rows = grepl("E\\. (ferrisi|falciformis)", Analysis)
        )
    )

# Save using your beautiful function
save_table_all_formats(supp_table_S2, "supplementary_table_validation_stats_rf_lab")
cat("✅ Supplementary Table S2 saved in all formats!\n")
