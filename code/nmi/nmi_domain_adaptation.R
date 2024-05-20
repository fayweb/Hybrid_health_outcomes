# Load necessary libraries
library("gtsummary")
library("glmnet")
library("dplyr")
library("ggplot2")

# Prepare data: Selecting necessary columns
genes <- hm %>%
    ungroup() %>%
    select(origin, all_of(Genes_v))

# Summary table by origin
table_summary <- genes %>%
    tbl_summary(by = origin, digits = everything() ~ 2) %>%
    add_p() %>%
    add_q()
print(table_summary)

# Rename and recode domain variable for modeling
genes <- genes %>%
    rename(domain = origin) %>%
    mutate(domain_binary = as.numeric(as.factor(domain)) - 1)

# Logistic regression to estimate domain membership probabilities
x <- as.matrix(genes[, -which(names(genes) == "domain")])
y <- genes$domain_binary
fit <- glmnet(x, y, family = 'binomial')

# Cross-validation to find optimal lambda
cv_fit <- cv.glmnet(x, y, family = 'binomial')

# Extract coefficients at the best lambda
optimal_lambda <- cv_fit$lambda.min
coefficients <- coef(fit, s = "lambda.min")

# Print coefficients to check
print(coefficients)

# Filter and prepare data for non-zero coefficients plotting
non_zero_coef <- coefficients[coefficients != 0]
data_to_plot <- data.frame(Feature = names(non_zero_coef), Importance = as.numeric(non_zero_coef))

# Plot feature importance
ggplot(data_to_plot, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flips the axes for easier reading
    labs(title = "Feature Importance", x = "Features", y = "Importance")

# Prepare for prediction and validation
set.seed(123)  # for reproducibility
train_indices <- sample(1:nrow(genes), size = 0.8 * nrow(genes))
train_data <- genes[train_indices, ]
test_data <- genes[-train_indices, ]

# Placeholder for response variable, replace with actual data
train_response <- rnorm(length(train_indices))

# Train model with weighted responses
train_fit <- glmnet(as.matrix(train_data[, -which(names(train_data) == "domain")]), 
                    train_response, weights = rep(1, length(train_response)), family = 'gaussian')

# Predict on test data using optimal lambda
test_predictions <- predict(train_fit, newx = as.matrix(test_data[, -which(names(test_data) == "domain")]), 
                            s = "lambda.min")

# Evaluate the model: Mean Squared Error
actual_responses <- rnorm(nrow(test_data))  # replace with actual test responses
mse <- mean((test_predictions - actual_responses)^2)
print(paste("Mean Squared Error on Test Data:", mse))

    