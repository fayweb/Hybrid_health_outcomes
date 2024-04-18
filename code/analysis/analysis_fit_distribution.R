# 7.3: We want to analyze the distribution of the predicted outcome variable
# "WL_max" (predicted weight loss dependent on the immune gene expression values)
# the distribution type is required for the downstream analysis
# It seems that our predicted weight loss variale fits a normal distribution


########## Analyzing the distribution of our data in order to 
# go on with the anaylsis 
Field %>% ggplot(aes(x = predicted_WL)) +
  geom_histogram(binwidth = 1.5)


##  predicted WL vs HI
Field %>%
    ggplot(aes(x = HI , y = predicted_WL , color = Sex)) +
    geom_smooth() +
    geom_point()

## body length vs predicted WL
Field %>%
    ggplot(aes(x = Body_Length , y = predicted_WL , color = Sex)) +
    geom_smooth() +
    geom_point()


## Let'S further analyse the distribution of WL
x <- Field$predicted_WL

descdist(data = x, discrete = FALSE)
descdist(data = x, discrete = FALSE, #data is continuous
         boot = 1000)

## 
normal_ <- fitdist(x, "norm")
weibull_ <- fitdist(x, "weibull")
gamma_ <- fitdist(x, "gamma")


## Now fit the distributions to the predicted weight loss data
tryDistrib(x, "normal") # -729.2451
tryDistrib(x, "binomial") # failed
tryDistrib(x, "student") # failed
tryDistrib(x, "weibull") # -740.7075
tryDistrib(x, "weibullshifted") # failed
# Seems like the normal distribution is a better fit in this case


## Compare again between normal and weibull
findGoodDist(x, "normal", "weibull") # normal here again


## plot the distributions
plot(normal_)
summary(normal_)
plot(gamma_)
summary(gamma_)
plot(weibull_)
summary(weibull_)


rm(normal_, weibull_, gamma_, )
