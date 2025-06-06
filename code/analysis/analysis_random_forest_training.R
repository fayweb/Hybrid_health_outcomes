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
saveRDS(WL_predict_gene, file =  paste0(cmodels, "WL_predict_gene.RDS"))
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
       dpi = 1000)


# Get variable importance from the WL_predict_gene fit
ImpData <- as.data.frame(randomForest::importance(WL_predict_gene))

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


