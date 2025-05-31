# ----------------------------------------------------------
# Script: Clean_after_impute.R
# Purpose: Final cosmetic cleaning for plots/legends
# Assumes all imputations, merges, and logical variables done earlier
# ----------------------------------------------------------

# Define infection group color mappings for ggplot
color_mapping <- c("Uninfected controls" = "#4ACAFF",
                   "E. ferrisi" = "#7A0092",
                   "E. falciformis" = "#FF0000")

color_mapping_f <- c("Uninfected" = "#4ACAFF",
                     "E. ferrisi" = "#7A0092",
                     "E. falciformis" = "#FF0000")

# Define ordered factor levels
factor_levels <- c("Uninfected controls", "E. ferrisi", "E. falciformis")
factor_levels_f <- c("Uninfected", "E. ferrisi", "E. falciformis")

# Clean and relabel parasite variables for plotting
hm$Parasite_primary <- gsub("_", ". ", hm$Parasite_primary)
hm$Parasite_challenge <- gsub("_", ". ", hm$Parasite_challenge)
hm$current_infection <- gsub("_", ". ", hm$current_infection)

hm$Parasite_primary <- gsub("uninfected", "Uninfected controls", hm$Parasite_primary)
hm$Parasite_challenge <- gsub("uninfected", "Uninfected controls", hm$Parasite_challenge)
hm$current_infection <- gsub("uninfected", "Uninfected controls", hm$current_infection)

hm$species_Eimeria <- gsub("_", ". ", hm$species_Eimeria)
hm$eimeriaSpecies <- gsub("_", ". ", hm$eimeriaSpecies)

hm$species_Eimeria <- gsub("uninfected", "Uninfected", hm$species_Eimeria)
hm$eimeriaSpecies <- gsub("uninfected", "Uninfected", hm$eimeriaSpecies)

# Apply factor levels
hm$Parasite_primary <- factor(hm$Parasite_primary, levels = factor_levels)
hm$Parasite_challenge <- factor(hm$Parasite_challenge, levels = factor_levels)
hm$current_infection <- factor(hm$current_infection, levels = factor_levels)

hm$eimeriaSpecies <- factor(hm$eimeriaSpecies, levels = factor_levels_f)
hm$species_Eimeria <- factor(hm$species_Eimeria, levels = factor_levels_f)

# Update plot-friendly labels
labels <- c("Uninfected controls", "*E. ferrisi*", "*E. falciformis*")
labels_f <- c("Uninfected", "*E. ferrisi*", "*E. falciformis*")

# Clean and order mouse strain factors for consistent plotting
hm$mouse_strain <- gsub("_", " ", hm$mouse_strain)
hm$mouse_strain <- factor(hm$mouse_strain, 
                          levels = names(sort(tapply(hm$WL_max, 
                                                     hm$mouse_strain, 
                                                     median))))

# Set lab subset for plotting convenience
lab <- hm %>% filter(origin == "Lab")

# Relevel for consistent legends
hm$immunization <- relevel(as.factor(hm$immunization), ref = "Uninfected controls")
hm$infection <- factor(hm$infection, levels = c("primary", "challenge"))

# Done – this script is now purely cosmetic and safe

