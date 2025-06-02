# 7.6: Could we derive tolerance out of the predicted health impact and 
# infection intenstities for each mouse? 
# If so, can we test the hybridicity on this derives tolerance variable?
field_infected <- Field %>%
    filter(infection_status == "TRUE" & !is.na(delta_ct_cewe_MminusE))
tolerance <- lm(formula = predicted_weight_loss ~  delta_ct_cewe_MminusE,
                data = field_infected)

summary(tolerance)

estimates <- ggpredict(tolerance)

plot(estimates)



