# 7.5: Exploring the relationships between variables and the predicted weight 
# loss of field mice
# Is there a relationship with the infection intensities with Eimeria spp.?
# Is there a relationship between bmi and predicted weight loss?
# We find a weak relationship of bmi and predicted weight loss

# Can the predicted weight loss be predicted by infection intensities
# infected samples
Field_infected <- Field %>%
  drop_na(delta_ct_cewe_MminusE) %>%
    filter(MC.Eimeria == TRUE)

ggplot(data = Field_infected, aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
    geom_point() +
    stat_smooth(method= "lm") 

# Is the predicted weight loss correlated to the infection intensities?
cor.test(Field_infected$predicted_WL, Field_infected$delta_ct_cewe_MminusE,
         method=c("pearson", "kendall", "spearman"))
# correlation non significant

# can the infection intensity predict weight loss?
model_WL <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field_infected)

summary(model_WL)
confint(model_WL)


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

Field2 <- Field %>%
  drop_na(OPG)

ggplot(data = Field2, aes(x = OPG, y = predicted_WL)) +
    geom_point() +
    stat_smooth(method= "lm") +
    scale_x_log10()


cor.test(Field2$predicted_WL, Field2$OPG)
model <- lm(predicted_WL ~  OPG, data = Field)


summary(model)
confint(model)



## Let'S dive into it further
# what about a combination of the oocyst and infection intensity data with qpcr?
model <- lm(predicted_WL ~  OPG * delta_ct_cewe_MminusE, data = Field)
summary(model)
confint(model)


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
    )

bmi <- lm(predicted_WL ~ BMI, data = Field)

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
cor.test(Field$BMI, Field$predicted_WL, use = "complete.obs")

# Summarize the linear model
summary(bmi)
# Confidence intervals for the model parameters
confint(bmi)


ggsave(filename = paste0(an_fi, "/BMI_WL.jpeg"),
       plot = bmi_plot, width = 5, height = 4, dpi = 1000)

