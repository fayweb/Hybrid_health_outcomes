

##All
fitWL_Sex <- parasiteLoad::analyse(data = Field,
                                   response = "predicted_WL",
                                   model = "normal",
                                   group = "Sex")

# Define or ensure MeanLoad is correctly defined. This is a placeholder.
MeanLoad <- function(L1, alpha, HI) {
    # Placeholder formula: adjust based on actual functional form
    return(L1 + alpha * HI)  # Adjust this formula as necessary
}

# Calculate predicted mean using MeanLoad for H0
# Assuming HI (host index or similar) is part of your data
predicted_mean_H0 <- MeanLoad(L1 = 9.3545283, alpha = -0.2786233, HI = Field$HI)  # Replace Field$HI with your actual data column

# Assuming normal distribution for simplicity; replace with actual distribution if different
predicted_values_H0 <- rnorm(length(predicted_mean_H0), mean = predicted_mean_H0, sd = 2.4458989)

# Calculate residuals
observed_values <- Field$predicted_WL  # Make sure this is your actual response variable
residuals_H0 <- observed_values - predicted_values_H0

# Histogram of residuals
hist(residuals_H0, main="Histogram of Residuals for H0", xlab="Residuals")

# Shapiro-Wilk test for normality
shapiro.test(residuals_H0)


# Plotting residuals against fitted values
plot(predicted_values_H0, residuals_H0, main="Residual vs. Fitted Plot for H0",
     xlab="Fitted Values", ylab="Residuals")
abline(h=0, col="red")  # Adding a horizontal line at 0

# Q-Q plot for normality
qqnorm(residuals_H0)
qqline(residuals_H0, col="red")
