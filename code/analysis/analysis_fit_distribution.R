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


# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
  # deals with fitdistr error:
  fit <- 
      tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
  return(list(fit = fit,
              loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"), 
              AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}


findGoodDist <- function(x, distribs, distribs2){
    l =lapply(distribs, function(i) tryDistrib(x, i))
    names(l) <- distribs
    print(l)
    listDistr <- lapply(distribs2, function(i){
        if (i %in% "t"){
            fitdistrplus::fitdist(x, i, start = list(df =2))
        } else {
            fitdistrplus::fitdist(x,i)
        }}
    ) 
    par(mfrow=c(2,2))
    denscomp(listDistr, legendtext=distribs2)
    cdfcomp(listDistr, legendtext=distribs2)
    qqcomp(listDistr, legendtext=distribs2)
    ppcomp(listDistr, legendtext=distribs2)
    par(mfrow=c(1,1))
}


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
