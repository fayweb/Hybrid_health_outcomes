# Load cleaned field data
field <- read.csv(paste0(dfield_inter, "/field_cleaned_intermediate.csv"))

# ---- Add infection intensities ----
ii <- read.csv("data/field/raw/CEWE_FECES_infection_intensities.txt")
ii$Mouse_ID <- gsub("AA_", "AA", ii$Mouse_ID)
field <- field %>% left_join(ii, by = intersect(colnames(field), colnames(ii)))
rm(ii)

# ---- Add amplicon species calls ----
amplicon <- read.csv("data/field/raw/Sample_selection_Metabarcoding_Complete.csv")
eimer_sp <- amplicon %>%
    dplyr::select(Mouse_ID, Species) %>%
    rename(amplicon_species = Species) %>%
    drop_na(amplicon_species)
eimer_sp$Mouse_ID <- gsub("_", "", eimer_sp$Mouse_ID)
field <- field %>% left_join(eimer_sp, by = "Mouse_ID")

# ---- Create species_Eimeria ----
field <- field %>%
    mutate(species_Eimeria = case_when(
        is.na(eimeriaSpecies) ~ amplicon_species,
        !is.na(eimeriaSpecies) ~ eimeriaSpecies
    ))
field$species_Eimeria <- gsub("Negative", "uninfected", field$species_Eimeria)
field$species_Eimeria <- as.factor(field$species_Eimeria)

# ---- Create infection_status ----
field$MC.Eimeria <- as.factor(field$MC.Eimeria)
field <- field %>%
    mutate(infection_status = case_when(
        is.na(MC.Eimeria) & species_Eimeria == "uninfected" ~ "FALSE",
        is.na(MC.Eimeria) & species_Eimeria %in% c("E_falciformis", "E_ferrisi") ~ "TRUE",
        !is.na(MC.Eimeria) ~ as.character(MC.Eimeria)
    ))
field$infection_status <- as.factor(field$infection_status)

rm(amplicon, eimer_sp)

# ---- Save final cleaned output ----
write.csv(field, paste0(dfield_final, "/field_cleaned_data.csv"), row.names = FALSE)

