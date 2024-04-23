
# add missing infection intensities
ii <- read.csv("data/field/raw/CEWE_FECES_infection_intensities.txt")

ii$Mouse_ID <- gsub(pattern = "AA_", replacement = "AA", x = ii$Mouse_ID)



field <- field %>%
    left_join(ii, by = intersect(colnames(field), colnames(ii)))


rm(ii)


amplicon <- read.csv("data/field/raw/Sample_selection_Metabarcoding_Complete.csv")



eimer_sp <- amplicon %>%
    dplyr::select(c(Mouse_ID, Species)) %>%
    rename(amplicon_species = Species) %>%
    drop_na(amplicon_species)

eimer_sp$Mouse_ID <- gsub(pattern = "_", replacement = "", x = eimer_sp$Mouse_ID)


field <- field %>%
    left_join(eimer_sp, by = "Mouse_ID")


# create new variable depending on amplicon sequencing
field <- field %>%
    mutate(species_Eimeria = case_when(
        is.na(eimeriaSpecies) ~ amplicon_species,
        !is.na(eimeriaSpecies) ~ eimeriaSpecies))

field$species_Eimeria <- gsub(pattern = "Negative", replacement = "uninfected", 
                              x = field$species_Eimeria)

field$species_Eimeria <- as.factor(field$species_Eimeria)
field$MC.Eimeria <- as.factor(field$MC.Eimeria)

field <- field %>%
    mutate(infection_status = case_when(
        is.na(MC.Eimeria) & species_Eimeria == "uninfected" ~ "FALSE",
        is.na(MC.Eimeria) & species_Eimeria == "E_falciformis" ~ "TRUE",
        is.na(MC.Eimeria) & species_Eimeria == "E_ferrisi" ~ "TRUE",
        !is.na(MC.Eimeria) ~ MC.Eimeria
    ))

field$infection_status <- as.factor(field$infection_status)

rm(amplicon, eimer_sp)
