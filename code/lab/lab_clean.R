## lab
## Change the strains of the parasite to the parasite strain names

lab <- lab %>%
    dplyr::mutate(Parasite_primary = case_when(
        primary_infection == "E64" ~ "E_ferrisi",
        primary_infection == "E88" ~ "E_falciformis",
        primary_infection == "Eflab" ~ "E_falciformis",
        primary_infection == "E139" ~ "E_ferrisi",
        primary_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))

lab <- lab %>%
    dplyr::mutate(Parasite_challenge = case_when(    
        challenge_infection == "E64" ~ "E_ferrisi",
        challenge_infection == "E88" ~ "E_falciformis",
        challenge_infection == "Eflab" ~ "E_falciformis",
        challenge_infection == "E139" ~ "E_ferrisi",
        challenge_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))

lab <- lab %>%
    dplyr::mutate(infection_history = case_when(
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "uninfected" ~ "uninfected",
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "E_ferrisi" ~ "uninfected_ferrisi",
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "E_falciformis" ~ "uninfected_falciformis",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "E_falciformis" ~ "falciformis_falciformis",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "E_ferrisi" ~ "falciformis_ferrisi",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "uninfected" ~ "falciformis_uninfected",
        Parasite_primary == "E_ferrisi" & 
            Parasite_challenge == "E_falciformis" ~ "ferrisi_falciformis",
        Parasite_primary == "E_ferrisi" & 
            Parasite_challenge == "E_ferrisi" ~ "ferrisi_ferrisi",
        Parasite_primary == "E_ferrisi" &
            Parasite_challenge == "uninfected" ~ "ferrisi_uninfected",
        TRUE ~ ""))


lab[sapply(lab, is.infinite)] <- NA


# IF a mouse dies before the end of the experiment the weight on the dpis after
# death will be NA. For further analysis I am changing those dpis to NA as well
# so that I can get the max dpi for each mouse

# replace the 0 with NA (incorrect experimental measurement)
lab$weight[lab$weight == 0] <- NA
lab$relative_weight[lab$relative_weight == 0] <- NA

# change dpis to NA after death
lab$dpi[is.na(lab$weight)] <- NA

### Add the variable end weight (relative weight at day of sacrifice)
# start by adding the variable dpi_max which inficates the last day of each mouse
lab <- lab%>%
    dplyr::filter(!weight == "NA") %>% 
    dplyr::group_by(EH_ID, infection) %>%
    dplyr::mutate(max_dpi = max(dpi), WL_max = (100 - min(relative_weight)))



#There are two measuremts for CXCR3
# We want to here keep the CXCR3_bio 
lab <- lab %>% 
    dplyr::select(-CXCR3)

#Now rename the CXCR3_bio to CXCR3
lab <- lab %>%
    dplyr::rename(CXCR3 = CXCR3_bio)


lab <- lab %>%
    dplyr::mutate(origin = "Lab")


# Rename column names to match each other
lab <- lab %>% 
    dplyr::rename(Mouse_ID = EH_ID, delta_ct_cewe_MminusE = delta, 
                  MC.Eimeria = Eim_MC, Feces_Weight = feces_weight)


# Here I create a new column, where we get the actual infection status
# According to the melting curve for eimeria 
lab <- lab %>%
    dplyr::mutate(current_infection = case_when(
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "TRUE" ~ "E_ferrisi",
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "uninfected" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "uninfected" & MC.Eimeria == "FALSE" ~ "uninfected",
        TRUE ~ Parasite_challenge
    ))


lab <- lab %>%
    dplyr::mutate(immunization = case_when(
        infection_history == "falciformis_ferrisi" ~ "heterologous",
        infection_history == "ferrisi_falciformis" ~ "heterologous",
        infection_history == "falciformis_uninfected" ~ "uninfected",
        infection_history == "ferrisi_uninfected" ~ "uninfected",
        infection_history == "ferrisi_ferrisi" ~ "homologous",
        infection_history == "falciformis_falciformis" ~ "homologous",
        infection_history == "uninfected_falciformis" ~ "naive",
        infection_history == "uninfected_ferrisi" ~ "naive",
        infection_history == "uninfected" ~ "uninfected",
        TRUE ~ "NA"
    ))

# Here I create a new column, where we get the actual infection status
# According to the melting curve for eimeria 
lab <- lab %>%
    dplyr::mutate(current_infection = case_when(
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "TRUE" ~ "E_ferrisi",
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "uninfected" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "uninfected" & MC.Eimeria == "FALSE" ~ "uninfected",
        TRUE ~ Parasite_challenge
    ))

lab <- lab %>%
    dplyr::mutate(immunization = case_when(
        infection_history == "falciformis_ferrisi" ~ "heterologous",
        infection_history == "ferrisi_falciformis" ~ "heterologous",
        infection_history == "falciformis_uninfected" ~ "uninfected",
        infection_history == "ferrisi_uninfected" ~ "uninfected",
        infection_history == "ferrisi_ferrisi" ~ "homologous",
        infection_history == "falciformis_falciformis" ~ "homologous",
        infection_history == "uninfected_falciformis" ~ "naive",
        infection_history == "uninfected_ferrisi" ~ "naive",
        infection_history == "uninfected" ~ "uninfected",
        TRUE ~ "NA"
    ))

Challenge <- lab


#### create output file:
write.csv(Challenge, paste0(dlab_final, "/lab_cleaned_data.csv"), row.names = FALSE)

