

################## hybrid effect
Field <- Field %>%
    mutate(HE = 2*HI*(1-HI)) #linearize HI) 
# tolerance = health impact / infection intensity

i <- Field %>%
    filter(Sex == "M") %>%
    drop_na(tolerance)

cor.test(i$HI, i$tolerance)
cor.test(Field$HI, Field$predicted_WL)

i <- Field %>%
    filter(Sex == "F") %>%
    drop_na(tolerance)

cor.test(i$HI, i$tolerance)
cor.test(i$HI, i$predicted_WL)


ggplot(Field, aes(x = HI, HE)) +
    geom_point() +
    geom_line()

ggplot(Field, aes(x = HE, predicted_WL, color = Sex)) +
    geom_smooth(method = lm, se = TRUE) 

ggplot(Field, aes(x = HE, tolerance, color = Sex)) +
    geom_jitter() +
    geom_smooth(method = lm, se = TRUE) 


ggplot(Field, aes(x = HE, OPG)) +
    geom_point() +
    geom_line()

model <- lm(formula = predicted_WL ~  HI + HE + HI * HE * Sex, data = Field)
summary(model)

df <- Field %>%
    filter(MC.Eimeria == TRUE)
###########################################################
#########################################################
model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

df %>%
    ggplot(aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
    geom_abline()

cor.test(df$predicted_WL, df$delta_ct_cewe_MminusE)

model <- lm(tolerance ~ HI + HE * Sex, df)
summary(model)


######################################################################
Field_asp <- Field %>%
    drop_na(Aspiculuris_sp) %>%
    filter(Aspiculuris_sp != 0) %>%
    mutate(tolerance = predicted_WL / Aspiculuris_sp) %>%
    drop_na(tolerance)

model <- lm(tolerance ~ Aspiculuris_sp, data = Field_asp)
summary(model)



model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

model <- lm(tolerance ~ HI + HE * Sex, df)
summary(model)
####################################################################
Field_asp <- Field %>%
    drop_na(Aspiculuris_sp) %>%
    filter(Aspiculuris_sp != 0) %>%
    mutate(tolerance = predicted_WL / Aspiculuris_sp) %>%
    drop_na(tolerance)

model <- lm(tolerance ~ Aspiculuris_sp, data = Field_asp)
summary(model)

#########################################################################
df <- Field %>%
    filter(eimeriaSpecies == "E_falciformis")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)


#################################
df <- Field %>%
    filter(eimeriaSpecies == "E_ferrisi")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

#################
df <- Field %>%
    filter(eimeriaSpecies == "E_falciformis")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)


#################################
df <- Field %>%
    drop_na(ILWE_Crypto_Ct) %>%
    filter(ILWE_Crypto_Ct != 0) %>%
    mutate(tolerance = predicted_WL / -ILWE_Crypto_Ct) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ ILWE_Crypto_Ct, 
                      data = df)

summary(model_tolerance)

###################################
df <- Field %>%
    drop_na(Syphacia_sp) %>%
    filter(Syphacia_sp != 0) %>%
    mutate(tolerance = predicted_WL / Syphacia_sp) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Syphacia_sp, 
                      data = df)

summary(model_tolerance)

#####################################
###################################
df <- Field %>%
    drop_na(Heterakis_sp) %>%
    filter(Heterakis_sp != 0) %>%
    mutate(tolerance = predicted_WL / Heterakis_sp) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Heterakis_sp, 
                      data = df)

summary(model_tolerance)


################################
df <- Field %>%
    drop_na(Trichuris_muris) %>%
    filter(Trichuris_muris != 0) %>%
    mutate(tolerance = predicted_WL / Trichuris_muris) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Trichuris_muris, 
                      data = df)

summary(model_tolerance)

##################################################################
##
Field <- Field %>%
    mutate(infection_intensity_Eim = delta_ct_cewe_MminusE)
Field_par <- Field %>%
    mutate(
        infected_Aspiculuris = Aspiculuris_sp != 0,
        infected_syphasia = Syphacia_sp != 0,
        infected_crypto = ILWE_Crypto_Ct != 0
    )

Field_par <- Field_par %>%
    mutate(
        infected_Aspiculuris = as.factor(infected_Aspiculuris),
        infected_syphasia = as.factor(infected_syphasia),
        infected_crypto = as.factor(infected_crypto)
    )

modelA <- lm(predicted_WL ~ MC.Eimeria + infected_Aspiculuris + HI + HE
             + infected_syphasia + infected_crypto, data = Field_par)

summary(modelA)

plot_summs(modelA, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue") -> plot1
plot1

ggsave(filename = "figures/coefficient_plot_model1.jpeg", plot = plot1, 
       width = 5, height = 4, dpi = 300)

#######################################
modelB <- lm(predicted_WL ~ MC.Eimeria * HE, Field)
summary(modelB)
plot_summs(modelB, plot.distributions = TRUE) 


## hybrid index + infectopm
Field$MC.Eimeria <- as.factor(Field$MC.Eimeria)
modelC <- lm(predicted_WL ~ MC.Eimeria *infection_intensity_Eim * HE + 
                 HI + HE, Field)
summary(modelC)

plot_summs(modelC, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue") -> plot2

plot2 

ggsave(filename = "figures/coefficient_plot_model2.jpeg", plot = plot2, 
       width = 5, height = 4, dpi = 300)

ggpredict(modelC, terms = c("MC.Eimeria", "infection_intensity_Eim")) %>%
    plot()


## hybrid index + infectopm
modelD <- lm(predicted_WL ~  HI + HE, Field)
summary(modelD)
plot_summs(modelD,  plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue")  -> plot_4

ggsave(filename = "figures/coefficient_plot_HE.jpeg", plot = plot_4, 
       width = 5, height = 4, dpi = 300)

plot_summs(modelC, modelD,  robust = TRUE, 
           scale = TRUE) -> model1_2
model1_2

##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
stargazer(modelA, modelB, modelC, modelD,
          type = "text",
          out = "tables/", 
          title = "Linear models - Contributions to predicted maximum weight loss",
          align = TRUE)stargazer_field_modelaA_D.doc

ggsave(filename = "figures/coefficient_plot_model1_2.jpeg", plot = model1_2, 
       width = 8, height = 6, dpi = 300)

model4 <- lm(predicted_WL ~  MC.Eimeria, Field)
summary(model4)

export_summs(modelA, modelB, modelC, modelD, scale = TRUE, to.file = "docx", 
             file.name = "tables/field_modelaA_D.docx")

summary(modelA)
summary(modelB)
summary(modelC)
summary(modelD)

plot_coefs(modelA, modelB, modelC, modelD, plot.distributions = TRUE)


##################################
#raincloud plots
ggplot(Field, aes(x = MC.Eimeria, y = predicted_WL, color = MC.Eimeria)) +
    geom_violin(trim = FALSE, alpha = 0.5) + 
    geom_quasirandom(aes(color = MC.Eimeria), size = 1, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(title = "Raincloud Plot of predicted_WL by MC.Eimeria Group",
         x = "MC.Eimeria Group",
         y = "Predicted WL")

# Filter out NA values in MC.Eimeria
Field_filtered <- Field %>% filter(!is.na(MC.Eimeria))


model <- lm(predicted_WL ~ MC.Eimeria, Field)
summary(model)

# Define colors# Define colorsField
colors <- c("infected" = "purple", "uninfected" = "cornflowerblue")

Field_filtered <- Field_filtered %>%
    mutate(MC.Eimeria = if_else(MC.Eimeria == "TRUE", "infected", "uninfected"))

ggplot(Field_filtered, aes(y = MC.Eimeria, x = predicted_WL, fill = MC.Eimeria)) + 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) + 
    scale_fill_manual(values = colors) +
    geom_boxplot(
        width = .15, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    stat_dots(
        # ploting on left side
        side = "left",
        # adjusting position
        justification = 1.1,
        # adjust grouping (binning) of observations
        binwidth = 0.25,
        alpha = 0.5) +
    geom_point(
        shape = 95,
        size = 15,
        alpha = .2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    coord_cartesian(ylim = c(1.2, 2.9), clip = "off") +
    theme_minimal() +
    labs(y = "Infection status with Eimeria spp.", 
         x = "Predicted maximum weight loss") +
    theme(legend.position = "none") -> raincloud_plots__eimeria

raincloud_plots__eimeria

ggsave(plot = raincloud_plots__eimeria, filename = "figures/raincloud_eimeria.jpeg", 
       width = 6, 
       height = 5, dpi = 1000)

h <- (raincloud_plots__eimeria | coe)

###############################################
########################
# raincloud plots parasites
model_par <- lm(predicted_WL ~ ILWE_Crypto_Ct + Aspiculuris_sp + Syphacia_sp, data = Field)
summary(model_par)

Field_par <- Field_par %>%
    pivot_longer(cols = c("infected_Aspiculuris", "infected_syphasia", 
                          "infected_crypto"), names_to = "parasite",
                 values_to = "infection_status")  %>% 
    drop_na(infection_status)


ggplot(Field_par, aes(y = infection_status, x = predicted_WL, 
                      fill = infection_status)) + 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.7,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) + 
    scale_fill_manual(values = colors) +
    geom_boxplot(
        width = .15, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    geom_point(
        shape = 95,
        size = 15,
        alpha = .2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    stat_dots(
        # ploting on left side
        side = "left",
        # adjusting position
        justification = 1.1,
        # adjust grouping (binning) of observations
        binwidth = 0.25) +
    coord_cartesian(ylim = c(1.2, 2.9), clip = "off") +
    theme_minimal() +
    labs(y = "Infection status with Eimerai spp.", 
         x = "Predicted detrimental immune signature") +
    facet_wrap(~parasite) -> parasites

parasites

ggsave(plot = parasites, filename = "figures/raincloud_parasites.jpeg", 
       width = 10, 
       height = 4, dpi = 1000)

###############################################################
################################################

model <- lm(predicted_WL ~ eimeriaSpecies, data = Field)
summary(model)



amplicon <- read.csv("https://raw.githubusercontent.com/ferreira-scm/Eimeria_AmpSeq/master/data/Wild/Sample_selection_Metabarcoding_Complete.csv")


glimpse(amplicon)
glimpse(Field)


eimer_sp <- amplicon %>%
    dplyr::select(c(Mouse_ID, Species)) %>%
    rename(amplicon_species = Species) %>%
    drop_na(amplicon_species)

eimer_sp$Mouse_ID <- gsub(pattern = "_", replacement = "", x = eimer_sp$Mouse_ID)


Field <- Field %>%
    left_join(eimer_sp, by = "Mouse_ID")




# create new variable depending on amplicon sequencing
Field <- Field %>%
    mutate(species_Eimeria = case_when(
        is.na(eimeriaSpecies) ~ amplicon_species,
        !is.na(eimeriaSpecies) ~ eimeriaSpecies))

Field$species_Eimeria <- gsub(pattern = "Negative", replacement = "uninfected", 
                              x = Field$species_Eimeria)

Field <- Field %>%
    mutate(infection_status = case_when(
        is.na(MC.Eimeria) & species_Eimeria == "uninfected" ~ "FALSE",
        is.na(MC.Eimeria) & species_Eimeria == "E_falciformis" ~ "TRUE",
        is.na(MC.Eimeria) & species_Eimeria == "E_ferrisi" ~ "TRUE",
        !is.na(MC.Eimeria) ~ MC.Eimeria
    ))

Field$infection_status <- as.factor(Field$infection_status)
############

infected <- parasiteLoad::analyse(data = Field,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "infection_status")



plot_infected <- 
    bananaPlot(mod = infected$H3,
               data = Field,
               response = "predicted_WL",
               group = "infection_status") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 

plot_infected

##############
model <- lm(predicted_WL ~ species_Eimeria, data = Field)
summary(model)


ggpredict(model = model) %>%
    plot()



##All
only_inf <- Field %>%
    dplyr::filter(species_Eimeria == c("E_falciformis", "E_ferrisi"))

only_inf$species_Eimeria <- as.factor(only_inf$species_Eimeria)

fitWL_species_Diff <- parasiteLoad::analyse(data = only_inf,
                                            response = "predicted_WL",
                                            model = "normal",
                                            group = "species_Eimeria")



plot_species <- 
    bananaPlot(mod = fitWL_species_Diff$H3,
               data = only_inf,
               response = "predicted_WL",
               group = "species_Eimeria") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 


plot_species

A <- Field %>%
    dplyr::select(c(Mouse_ID, infection_status, species_Eimeria))


write.csv(A, "tables/table_for_saskia.csv", row.names = FALSE)

