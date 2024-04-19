# 7.6: Could we derive tolerance out of the predicted health impact and 
# infection intenstities for each mouse? 
# If so, can we test the hybridicity on this derives tolerance variable?
Field_tol <- Field %>%
    mutate(tolerance = predicted_WL / delta_ct_cewe_MminusE)

Field_tol <- Field_tol %>%
    filter(!is.na(tolerance), MC.Eimeria == TRUE)

summary(Field_tol$tolerance)

Field_tol <- Field_tol %>%
    filter(tolerance > -5, tolerance < 30)


summary(Field_tol$tolerance)

hist(Field_tol$tolerance)

Field_tol %>%
    ggplot(aes(tolerance)) +
    geom_histogram()

parasiteLoad::getParamBounds("normal", data = Field_tol, response = "tolerance")

x <- Field_tol$tolerance

tryDistrib(x, "normal") #-206.8664
tryDistrib(x, "binomial")
tryDistrib(x, "student")
tryDistrib(x, "weibull")
tryDistrib(x, "weibullshifted")

Field_tol$Sex <- as.factor(Field_tol$Sex)


##All
fitWL_tol <- parasiteLoad::analyse(data = Field_tol,
                                   response = "tolerance",
                                   model = "normal",
                                   group = "Sex")


##plot_tolerance_Sex <- 
#    bananaPlot(mod = fitWL_tol$H3,
 #              data = Field_tol,
  #             response = "tolerance",
   #            group = "Sex") +
    #scale_fill_manual(values = c("blueviolet", "limegreen")) +
    #scale_color_manual(values = c("blueviolet", "limegreen")) +
    #theme_bw() 


#plot_tolerance_Sex

# Create HI bar
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)),
                        aes(x=hi, y=1, fill = hi)) +
    geom_tile() +
    theme_void() +
    scale_fill_gradient(low = "blue", high = "red")  + 
    scale_x_continuous(expand=c(.01,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position = 'none')

