# Load intermediate
field <- read.csv(paste0(dfield_inter, "/field_imported_raw.csv"))

# Add origin flag
field <- field %>%
    mutate(origin = "Field")

# Adjust parasite names
field <- field %>%
    mutate(eimeriaSpecies = case_when(
        eimeriaSpecies == "Negative" ~ "uninfected",
        eimeriaSpecies == "" ~ "NA",
        TRUE ~ eimeriaSpecies
    ))

# Clean Mouse_ID formatting
field$Mouse_ID <- gsub("_", "", field$Mouse_ID)

# Save cleaned intermediate
write.csv(field, paste0(dfield_inter, "/field_cleaned_intermediate.csv"), row.names = FALSE)
