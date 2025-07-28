# ***********************************************************
# Part 7: Analysis                           ----
# ***********************************************************
#----------------------------------------------------------*
#7: Random forest 
# Training and testing a random forest that predicts weight los
# on experimental infections with Eimeria spp. 
# Requires: hm
# Creates random forest model: WL_predict_gene.RData
#----------------------------------------------------------*

lab <- hm %>%
    filter(origin == "Lab")
#select the imputed gene columns
gene_m <-  lab %>%
    ungroup() %>%
    dplyr::select(c(Mouse_ID, all_of(Genes_v), WL_max))

# select only the genes
genes <- gene_m %>%
  dplyr::select(-Mouse_ID)

# select the genes and the weight loss
gene_W <- lab  %>%
    ungroup() %>%
    dplyr::select(c(all_of(Genes_v), WL_max))

repeat_cv <- trainControl(method = "repeatedcv", #repeated cross validation
                           number = 5, # 5 fold cross validation
                           repeats = 3)

# split data into training and test
set.seed(333) # this will help us reproduce this random assignment

# in this way we can pick the random numbers
training.samples <- createDataPartition(y = gene_W$WL_max, p = .7, list = FALSE) 

# this is the partiicition! In this case 0.7 = training data and 0.3 = testing
# we don't want to get a list in return
train.data <- gene_W[training.samples, ] 
test.data <- gene_W[-training.samples, ] 


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

predict_WL_cv$fit.var.exp
predict_WL_cv$fit.mse
par(mar=c(1,1,1,1))

##################
##################
########## Plots

root_mean <- plot(predict_WL_cv)

# Root Mean Squared Error (observed vs. predicted) from each Bootstrap 
# iteration (cross-validation)
mean_error <- plot(predict_WL_cv, stat = "mse")

#Percent variance explained from specified fit model
model_var <- plot(predict_WL_cv, stat = "var.exp")

#Mean Absolute Error from each Bootstrapped model
abs_error <- plot(predict_WL_cv, stat = "mae")


#d# ---------------------------------------------------------------------------------------------------
error_random  <- plot(WL_predict_gene)

## ---------------------------------------------------------------------------------------------------
# number of trees with lowest MSE
which.min(WL_predict_gene$mse)

# RMSE of this optimal random forest
sqrt(WL_predict_gene$mse[which.min(WL_predict_gene$mse)])

WL_predict_gene$mtry
oob_error_rate <- WL_predict_gene$mse[WL_predict_gene$ntree]
oob_error_rate <- 1 - sum(diag(WL_predict_gene$confusion)) / sum(WL_predict_gene$confusion)


### Visualize variable importance ---
#Call importance() function on the model model to check how the attributes used 
# as predictors affect our WL_predict_gene
ImpData <- as.data.frame(randomForest::importance(WL_predict_gene))
ImpData$Var.Names <- row.names(ImpData)
varImp(WL_predict_gene)
var_imp <- as.data.frame(varImp(WL_predict_gene))
var_imp$Genes <- row.names(var_imp)
var_imp <- var_imp %>%
    rename(Importance = Overall)


# Assuming var_imp is your data frame with variables 'Importance' and 'Genes'
var_imp <- var_imp %>%
    mutate(Genes = factor(Genes, levels = Genes[order(-Importance)])) # Reorder Genes by decreasing Importance

# Create the plot with a color scale
 importance_plot <- 
     ggplot(var_imp, aes(x = reorder(Genes, Importance), y = Importance, fill = Importance)) +
        geom_col() + # Use geom_col for a bar plot; it's more appropriate for importance scores
        coord_flip() + # Flip the coordinates to make it easier to read
        labs(x = "Genes", y = "IncNodePurity") +#, title = "Variable Importance of Genes") +
        theme_minimal() + # A clean, minimal theme
        theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjust text angle for x-axis labels if needed
              legend.title = element_blank()) + # Remove the legend title if desired
        scale_fill_viridis_c(option = "magma", direction = -1) + # Apply a Viridis color scale with the 'magma' option
        theme(legend.position = c(0.8, 0.4))
 
 importance_plot
## S3 method for class 'randomForest'
plot(WL_predict_gene, type = "l", main=deparse(substitute(x)))

variable_importance <- varImpPlot(WL_predict_gene)


ggsave(filename = paste0(an_fi, "/variable_imporance_random.pdf"),
       #width = 6, height = 5, 
       dpi = 300)

ggsave(filename = paste0(an_fi, "/variable_imporance_random.jpeg"),
       #width = 6, height = 5, 
       dpi = 300)


# Get variable importance from the WL_predict_gene fit
ImpData <- as.data.frame(randomForest::importance(WL_predict_gene))



# ***********************************************************
# PARTIAL DEPENDENCE PLOTS - INSERT AFTER VARIABLE IMPORTANCE SAVING
# ***********************************************************

# Install and load pdp package if needed
if (!require(pdp)) {
    install.packages("pdp")
    library(pdp)
}

cat("Creating partial dependence plots for top 3 predictors...\n")

# Get top 3 predictors based on importance
top_3_genes <- var_imp %>% 
    arrange(desc(Importance)) %>%
    head(3) %>% 
    pull(Genes)

cat("Top 3 predictors:", paste(top_3_genes, collapse = ", "), "\n")

# Generate partial dependence plots
pdp_plots <- list()
for(gene in top_3_genes) {
    cat("Processing", gene, "...\n")
    
    # Calculate partial dependence
    pd_data <- partial(WL_predict_gene, pred.var = gene, train = train.data)
    
    # Create individual plot
    pdp_plots[[gene]] <- pd_data %>%
        autoplot() + 
        labs(title = gene,  # Shorter title for better fit
             x = paste(gene, "Expression"), 
             y = "Predicted Weight Loss (%)") +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 9))
}

# Combine partial dependence plots
pdp_combined <- wrap_plots(pdp_plots, ncol = 3)

# Add overall title
pdp_final <- pdp_combined + 
    plot_annotation(
        title = "",
        theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
    )

# Save partial dependence plots
ggsave(filename = paste0(an_fi, "/partial_dependence_plots.pdf"), 
       plot = pdp_final, width = 15, height = 5, dpi = 300)

ggsave(filename = paste0(an_fi, "/partial_dependence_plots.jpeg"), 
       plot = pdp_final, width = 15, height = 5, dpi = 300)

cat("✅ Partial dependence plots saved!\n")

# Create individual PDP plots for supplementary material
for(i in 1:length(top_3_genes)) {
    gene <- top_3_genes[i]
    ggsave(filename = paste0(an_fi, "/pdp_individual_", gene, ".pdf"), 
           plot = pdp_plots[[gene]], width = 5, height = 4, dpi = 300)
}

cat("✅ Individual PDP plots saved for supplementary material!\n")
cat("=====================================\n")




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


# what is the correlation between predicted and actual data?
cor(result$WL_max, result$predictions, 
    method = c("pearson", "kendall", "spearman"))

cor(result$WL_max, result$predictions, 
    method = "pearson")

model <- lm(predictions ~ WL_max, data = test_lab)
summary(model)

ggpredict(model, terms = c("WL_max")) %>% 
    plot(colors = "blue") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss during infections") +
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

ggsave(filename = paste0(an_fi, "/correlation_pred_obs_random.pdf"),
       width = 6, height = 5, dpi = 300)


# ***********************************************************
# Random Forest Model Diagnostics and Comparison
# Add this to your existing random forest script
# ***********************************************************
# ***********************************************************
# Part 1: Model Diagnostics for Test Set Model (43% model)
# ***********************************************************

# 1.1: Residuals vs Fitted Plot
residuals_data <- data.frame(
    fitted = result$predictions,
    residuals = result$WL_max - result$predictions,
    observed = result$WL_max
)

# Residuals vs Fitted
residuals_fitted_plot <- ggplot(residuals_data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.3) +
    labs(x = "Fitted Values", 
         y = "Residuals",
         title = "Residuals vs Fitted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# 1.2: Q-Q Plot of Residuals
qq_plot <- ggplot(residuals_data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(x = "Theoretical Quantiles", 
         y = "Sample Quantiles",
         title = "Q-Q Plot of Residuals") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# 1.3: Cross-validation Performance Curves
# Extract CV results from your existing predict_WL_cv object
cv_results <- data.frame(
    iteration = 1:length(predict_WL_cv$y.rmse),
    rmse = predict_WL_cv$y.rmse,
    var_exp = predict_WL_cv$fit.var.exp,
    mae = predict_WL_cv$y.mae
)

# Quick fix - create diagnostic panel with the plots that worked
diagnostic_panel <- grid.arrange(
    residuals_fitted_plot, qq_plot,
    ncol = 2, nrow = 1,
    top = "Random Forest Model Diagnostics"
)


# RMSE across CV iterations
cv_rmse_plot <- ggplot(cv_results, aes(x = iteration, y = rmse)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = mean(cv_results$rmse), 
               linetype = "dashed", color = "red") +
    labs(x = "Cross-Validation Iteration", 
         y = "RMSE",
         title = "Cross-Validation RMSE") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Variance explained across CV iterations
cv_var_plot <- ggplot(cv_results, aes(x = iteration, y = var_exp)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = mean(cv_results$var_exp), 
               linetype = "dashed", color = "red") +
    labs(x = "Cross-Validation Iteration", 
         y = "Variance Explained (%)",
         title = "Cross-Validation Variance Explained") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Combine diagnostic plots
diagnostic_panel <- grid.arrange(
    residuals_fitted_plot, qq_plot,
    cv_rmse_plot, cv_var_plot,
    ncol = 2, nrow = 2
)

# Save diagnostic plots
ggsave(filename = paste0(an_fi, "/rf_model_diagnostics.pdf"),
       plot = diagnostic_panel, width = 12, height = 10, dpi = 300)
ggsave(filename = paste0(an_fi, "/rf_model_diagnostics.jpeg"),
       plot = diagnostic_panel, width = 12, height = 10, dpi = 300)

# ***********************************************************
# Part 2: Model Comparison Table
# ***********************************************************

# Create comparison table between test model and final model
model_comparison <- data.frame(
    Model = c("Test Set Model", "Final Model"),
    Dataset = c("Training subset (70%)", "Complete dataset"),
    n_samples = c(nrow(train.data), nrow(gene_W)),
    n_trees = c(26, 308),
    Variance_Explained = c(43.55, 47.6),
    MSE = c(43.23, 32.47),
    RMSE = c(sqrt(43.23), sqrt(32.47)),
    Test_Correlation = c(0.79, "Applied to wild mice"),
    Purpose = c("Validation & performance assessment", "Wild mouse prediction")
)

# Print the table
print("Random Forest Model Comparison:")
print(model_comparison)

# Save as formatted table
write.csv(model_comparison, 
          file = paste0(tables, "/rf_model_comparison.csv"), 
          row.names = FALSE)

# Create a formatted display table
comparison_table <- model_comparison %>%
    kbl(caption = "Random Forest Model Comparison",
        digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(1, background = "#f0f0f0") %>%
    row_spec(2, background = "#e6f3ff")

# Save formatted table
save_kable(comparison_table, 
           file = paste0(tables, "/rf_model_comparison_formatted.html"))

# ***********************************************************
# Part 3: Top Predictor Correlation Matrix
# ***********************************************************

# Get top 5 predictors based on importance
top_predictors <- var_imp %>%
    top_n(5, Importance) %>%
    pull(Genes)

# Create correlation matrix for top predictors
top_pred_data <- gene_W %>%
    dplyr::select(all_of(top_predictors))

# Calculate correlation matrix
cor_matrix <- cor(top_pred_data, use = "complete.obs")

# Create correlation plot with magma color scheme (matching your importance plot)
pdf(paste0(an_fi, "/top_predictors_correlation.pdf"), width = 8, height = 6)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         addCoef.col = "white",  # Changed to white for better contrast
         number.cex = 0.8,
         title = "Correlation Matrix: Top 5 Predictors",
         mar = c(0,0,2,0),
         col = viridis::magma(200, direction = -1))  # Add magma colors
dev.off()

# Alternative: Create a custom color palette that matches magma
magma_colors <- magma(200, direction = -1)

pdf(paste0(an_fi, "/top_predictors_correlation_magma.pdf"), width = 8, height = 6)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         addCoef.col = "white",
         number.cex = 0.8,
         title = "Correlation Matrix: Top 5 Predictors",
         mar = c(0,0,2,0),
         col = magma_colors)
dev.off()

# Save as JPEG with proper dimensions
jpeg(paste0(an_fi, "/top_predictors_correlation_magma.jpeg"), 
     width = 8 * 300, height = 6 * 300, res = 300)  # Need to specify width/height in pixels for jpeg
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         addCoef.col = "white",
         number.cex = 0.8,
         title = "Correlation Matrix: Top 5 Predictors",
         mar = c(0,0,2,0),
         col = magma_colors)
dev.off()

# Original ggplot version that worked
cor_data <- as.data.frame(as.table(cor_matrix))
names(cor_data) <- c("Gene1", "Gene2", "Correlation")

cor_heatmap <- ggplot(cor_data, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "", title = "Top 5 Predictors: Correlation Matrix") +
    geom_text(aes(label = round(Correlation, 2)), size = 3)

ggsave(filename = paste0(an_fi, "/top_predictors_correlation_ggplot.jpeg"),
       plot = cor_heatmap, width = 8, height = 6, dpi = 300)


### prdiction by group
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
       dpi = 300)

ggsave(plot = predictions_random_for_lab, 
       filename = paste0(an_fi, "/predictions_random_for_lab.jpeg"), 
       width = 6, height = 5,
       dpi = 300)
# ***********************************************************
# Part 4: Summary Statistics
# ***********************************************************

# Model performance summary
cat("\n=== RANDOM FOREST MODEL SUMMARY ===\n")
cat("Test Set Performance (n =", nrow(test.data), "):\n")
cat("- Correlation (r):", round(cor(result$WL_max, result$predictions), 3), "\n")
cat("- R-squared:", round(cor(result$WL_max, result$predictions)^2, 3), "\n")
cat("- RMSE:", round(sqrt(mean((result$WL_max - result$predictions)^2)), 3), "\n")

cat("\nCross-Validation Results:\n")
cat("- Mean RMSE:", round(mean(cv_results$rmse), 3), "\n")
cat("- Mean Variance Explained:", round(mean(cv_results$var_exp), 3), "%\n")

cat("\nFinal Model (n =", nrow(gene_W), "):\n")
cat("- Variance Explained:", WL_predict_gene$rsq[length(WL_predict_gene$rsq)] * 100, "%\n")
cat("- MSE:", round(WL_predict_gene$mse[length(WL_predict_gene$mse)], 3), "\n")

cat("\nTop 5 Predictors:\n")
top_5 <- var_imp %>% 
    arrange(desc(Importance)) %>%
    head(5)
for(i in 1:nrow(top_5)) {
    cat(paste0(i, ". ", top_5$Genes[i], " (", round(top_5$Importance[i], 1), ")\n"))
}

cat("=====================================\n")

# ***********************************************************
# Part 5: Model Justification Text
# ***********************************************************

# Generate text for methods/results
cat("\n=== TEXT FOR MANUSCRIPT ===\n")
cat("Model Development Strategy:\n")
cat("We employed a two-stage random forest approach. First, we split the laboratory dataset")
cat("into training (70%) and testing (30%) sets to assess predictive performance and avoid")
cat("overfitting. After confirming robust performance (r = 0.79), we retrained the model")
cat("on the complete laboratory dataset to maximize predictive power for application to")
cat("wild-caught mice. This approach balances validation rigor with practical utility.\n")
cat("============================\n")



########################################################################
## Decision to re-train the random forest to the complete data set before 
# application to wild samples
set.seed(232)

#train the model
WL_predict_gene <- randomForest(WL_max ~., data = gene_W, 
                                proximity = TRUE, ntree = 308) 

# ntree = number of trees     
# save the model 
saveRDS(WL_predict_gene, file =  paste0(cmodels, "WL_predict_gene.RDS"))
print(WL_predict_gene)

predict_WL_cv <- rf.crossValidation(x = WL_predict_gene, xdata = gene_W, 
                                    p = 0.10, n = 99, ntree = 308)

predict_WL_cv$fit.var.exp
predict_WL_cv$fit.mse
par(mar=c(1,1,1,1))

##################
##################
########## Plots

root_mean <- plot(predict_WL_cv)

# Root Mean Squared Error (observed vs. predicted) from each Bootstrap 
# iteration (cross-validation)
mean_error <- plot(predict_WL_cv, stat = "mse")

#Percent variance explained from specified fit model
model_var <- plot(predict_WL_cv, stat = "var.exp")

#Mean Absolute Error from each Bootstrapped model
abs_error <- plot(predict_WL_cv, stat = "mae")


#d# ---------------------------------------------------------------------------------------------------
error_random  <- plot(WL_predict_gene)

## ---------------------------------------------------------------------------------------------------
# number of trees with lowest MSE
which.min(WL_predict_gene$mse)

# RMSE of this optimal random forest
sqrt(WL_predict_gene$mse[which.min(WL_predict_gene$mse)])

WL_predict_gene$mtry
oob_error_rate <- WL_predict_gene$mse[WL_predict_gene$ntree]
oob_error_rate <- 1 - sum(diag(WL_predict_gene$confusion)) / sum(WL_predict_gene$confusion)
 

