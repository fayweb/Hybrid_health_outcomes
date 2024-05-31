# 7.5: Exploring the relationships between variables and the predicted weight 
# loss of field mice
# Is there a relationship with the infection intensities with Eimeria spp.?
# Is there a relationship between bmi and predicted weight loss?
# We find a weak relationship of bmi and predicted weight loss

# Can the predicted weight loss be predicted by infection intensities
# infected samples
Field_infected <- Field %>%
  drop_na(delta_ct_cewe_MminusE) %>%
    filter(infection_status == TRUE)

ggplot(data = Field_infected, aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
    geom_point() +
    stat_smooth(method= "lm") 


# Create the scatter plot
ggplot(data = Field_infected, aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
    geom_point(color = "black", size = 3, alpha = 0.6) +   # Plot the data points with specified color, size, and transparency
    geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", size = 1) +  # Add a linear trend line with confidence interval
    geom_smooth(method = "loess", se = TRUE, color = "red", linetype = "dashed", size = 1) +  # Add a smoothing line with confidence interval
    labs( x = "Infection Intensity of Eimeria spp.",
        y = "Predicted Maximum Weight Loss",) +
    theme_minimal(base_size = 15) +  # Use a minimal theme with larger base font size for better readability
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center the plot title and adjust its style
        plot.subtitle = element_text(hjust = 0.5, size = 14),  # Center the plot subtitle and adjust its size
        axis.title = element_text(face = "bold"),  # Bold axis titles
        panel.grid.major = element_line(color = "grey80"),  # Subtle grid lines for major grids
        panel.grid.minor = element_blank()  # Remove minor grid lines
    ) -> Cor_infection_wl

Cor_infection_wl

ggsave(paste0(an_fi, "/correlation_intensity_pred_WL.jpeg"), Cor_infection_wl)


plot(Field_infected$predicted_WL, Field_infected$delta_ct_cewe_MminusE)
shapiro.test(Field_infected$predicted_WL) # normal distribution
shapiro.test(Field_infected$delta_ct_cewe_MminusE) # non normal distribution

lm_model <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
               data = Field_infected)
summary(lm_model)
plot(lm_model$residuals)

# Is the predicted weight loss correlated to the infection intensities?
cor.test(Field_infected$predicted_WL, Field_infected$delta_ct_cewe_MminusE,
         method=c("spearman"), exact=FALSE)
# correlation non significant

# can the infection intensity predict weight loss?
model_WL <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field_infected)
model_fa <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field_infected %>%
                   filter(species_Eimeria == "E. falciformis"))
model_fe <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field_infected %>%
                   filter(species_Eimeria == "E. ferrisi"))
summary(model_WL)
summary(model_fa)
summary(model_fe)
confint(model_WL)


fa <- Field_infected %>%
    filter(species_Eimeria == "E. falciformis")

cor.test(fa$predicted_WL, fa$delta_ct_cewe_MminusE,
         method=c("spearman"), exact=FALSE)
fe <- Field_infected %>%
    filter(species_Eimeria == "E. ferrisi")
cor.test(fe$predicted_WL, fe$delta_ct_cewe_MminusE,
         method=c("spearman"), exact=FALSE)

model_WL_infection <- 
    ggpredict(model_WL) %>% 
    plot(colors = "darkorange2") +   # Use a refined shade of blue
    labs(title = NULL) +  
    xlab("Infection intensity of Eimeria spp. in the caecum") +
    ylab("Predicted values of weight loss") +
    theme_bw()

model_WL_infection


## -How about the oocysts?
ggplot(data = Field, aes(x = OPG, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") +
  scale_x_log10()

# Create the scatter plot for the variable OPG
Cor_OPG_wl <- ggplot(data = Field, aes(x = OPG, y = predicted_WL)) +
    geom_point(color = "black", size = 3, alpha = 0.6) +   # Plot the data points with specified color, size, and transparency
    geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", size = 1) +  # Add a linear trend line with confidence interval
    geom_smooth(method = "loess", se = TRUE, color = "red", linetype = "dashed", size = 1) +  # Add a smoothing line with confidence interval
    labs(
        x = "Predicted Maximum Weight Loss",
        y = "Infection Intensity of Eimeria spp. (OPG)"
    ) +
    theme_minimal(base_size = 15) +  # Use a minimal theme with larger base font size for better readability
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center the plot title and adjust its style
        plot.subtitle = element_text(hjust = 0.5, size = 14),  # Center the plot subtitle and adjust its size
        axis.title = element_text(face = "bold"),  # Bold axis titles
        panel.grid.major = element_line(color = "grey80"),  # Subtle grid lines for major grids
        panel.grid.minor = element_blank()  # Remove minor grid lines
    ) +
    scale_x_log10()

Cor_OPG_wl

shapiro.test(Field$OPG) # non normal distribution


cor.test(Field$predicted_WL, Field$OPG,  method= "spearman", exact = FALSE)

model <- lm(predicted_WL ~  OPG * species_Eimeria, data = Field)

summary(model)
confint(model)



## Let'S dive into it further
# what about a combination of the oocyst and infection intensity data with qpcr?
model <- lm(predicted_WL ~  infection_status * delta_ct_cewe_MminusE, data = Field)
summary(model)
confint(model)

# Plot the summary with enhanced aesthetics
plot_summs(model, 
           scale = TRUE, 
           inner_ci_level = 0.95, 
           outer_ci_level = 0.99) +
    theme_minimal() +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
    ) +
    labs(
        x = "Estimate",
        y = "Predictor"
    ) +
    scale_color_manual(values = c("blue", "red")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") -> int_MC

ggsave(paste0(an_fi, "/intensity_melting_curve.jpeg"), int_MC)

## LetÄ's create the BMI variable and see how it works
Field <- Field %>%
    dplyr::mutate(BMI = Body_Weight / (Body_Length / 100)^2)

# plotting the bmi vs predicted weight loss
ggplot(data = Field, aes(x = BMI, y = predicted_WL)) +
    geom_point(color = "blue", alpha = 0.6) +  # Adjust point color and transparency
    geom_smooth(method= "lm", color = "red") +  # Highlight the linear model
    labs(
        title = "Relationship between BMI and Predicted WL",
        y = "Predicted Weight Loss",
        x = "Body Mass Index (BMI)"
    ) +
    theme_minimal() +  # Use a minimalistic theme
    theme(
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.text = element_text(color = "black"),  # Ensure text is readable
        legend.position = "none"  # Hide the legend if unnecessary
    ) -> bmi_plot

bmi_plot

bmi <- lm(predicted_WL ~ BMI, data = Field)
shapiro.test(Field$BMI) # non normal distribution

# fitting the bmi model
predictions <- ggpredict(bmi, terms = c("BMI"))

# Use ggplot for plotting the predictions
ggplot(predictions, aes(x = x, y = predicted)) +
    geom_line(color = "blue") +  # Change to geom_line() for a continuous prediction
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +  # Confidence interval
    labs(
        title = "Predicted weight loss based on BMI",
        y = "Predicted Health Impact",
        x = "Body Mass Index (BMI)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Perform correlation test
cor.test(Field$BMI, Field$predicted_WL, method = "spearman", use = "complete.obs")

# Summarize the linear model
summary(bmi)
# Confidence intervals for the model parameters
confint(bmi)


ggsave(filename = paste0(an_fi, "/BMI_WL.jpeg"),
       plot = bmi_plot, width = 5, height = 4, dpi = 1000)

