# Load raw field data
field <- read.csv(paste0(dfield_raw, "/Field_infection_data.csv"))

# Save as intermediate (untouched)
write.csv(field, paste0(dfield_inter, "/field_imported_raw.csv"), row.names = FALSE)
