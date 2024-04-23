field <- field %>%
    dplyr::mutate(origin = "Field")

# Adjust the parasite names to fit the lab
field <- field %>%
    dplyr::mutate(eimeriaSpecies = case_when(
        eimeriaSpecies == "Negative" ~ "uninfected",
        eimeriaSpecies == "" ~ "NA",
        TRUE ~ eimeriaSpecies))

field$Mouse_ID <- gsub("_", "", field$Mouse_ID)

#### create output file:
write.csv(field, paste0(dfield_final, "/field_cleaned_data.csv"), row.names = FALSE)

