# 7.6: Could we derive tolerance out of the predicted health impact and 
# infection intenstities for each mouse? 
# If so, can we test the hybridicity on this derives tolerance variable?
tolerance <- lm(predicted_WL ~  delta_ct_cewe_MminusE, Field %>% 
                    filter(infection_status == "TRUE"))

summary(tolerance)

estimates <- ggpredict(tolerance)

plot(estimates)



