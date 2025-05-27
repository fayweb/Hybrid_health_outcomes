CXCL9_lab <- lm(WL_max ~ CXCL9, data = lab)
summary(CXCL9_lab)

CXCL9_field <- lm(predicted_WL ~ CXCL9, data = Field)
summary(CXCL9_field)

# combine "weight loss" 
# Add the weight loss colum to the field
Field <- Field %>% mutate(WL_max = predicted_WL)

hm_WL <- rbind(Field, lab)

CXCL9_combined <- lm(WL_max ~ CXCL9, data = hm_WL)
summary(CXCL9_combined)

# Train model on lab data only
lab_model <- lm(WL_max ~ CXCL9, data = lab)
summary(lab_model)

# Apply trained model to predict field data
field_predictions <- predict(lab_model, newdata = Field)

# Test how well lab-trained model predicts field outcomes
cor(field_predictions, Field$predicted_WL)

# Get correlation coefficient
correlation_result <- cor(field_predictions, Field$predicted_WL)
print(paste("Correlation:", correlation_result))

# For significance testing
cor.test(field_predictions, Field$predicted_WL)
