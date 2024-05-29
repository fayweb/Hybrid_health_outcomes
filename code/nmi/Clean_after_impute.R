# 4.5 Cleaning after imputing  
# Correct the parasite labels for legend
# Add custom colours for parasites throughout the scripts
# Creates factor levels for parasite strains
#----------------------------------------------------------*
# Use this to get the parasite colors in your ggplots:
#scale_color_manual(values = color_mapping) +

######### Decision: Removing IL.10 from gene selection due data largely missing
# remove il_10
#### vectors for selecting genes for analysis
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
               "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
               "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
               "TICAM1", "TNF")


# Creating colors and labels for following scripts 
# preparing the parasite names for figure legends
hm$Parasite_primary <- gsub(pattern = "_", ". ",  hm$Parasite_primary)

hm$Parasite_challenge <- gsub(pattern = "_", ". ",  hm$Parasite_challenge)

hm$Parasite_primary <- gsub(pattern = "uninfected", "Uninfected controls",  
                            hm$Parasite_primary)

hm$Parasite_challenge <- gsub(pattern = "uninfected", "Uninfected controls",  
                              hm$Parasite_challenge)

hm$current_infection <- gsub(pattern = "_", ". ",  hm$current_infection)

hm$current_infection <- gsub(pattern = "uninfected", "Uninfected controls",
                             hm$current_infection)

Challenge$Parasite_primary <- gsub(pattern = "_", ". ",  
                                   Challenge$Parasite_primary)

Challenge$Parasite_challenge <- gsub(pattern = "_", ". ",  
                                     Challenge$Parasite_challenge)

Challenge$Parasite_primary <- gsub(pattern = "uninfected", "Uninfected controls",  
                                   Challenge$Parasite_primary)

Challenge$Parasite_challenge <- gsub(pattern = "uninfected", "Uninfected controls",  
                                     Challenge$Parasite_challenge)

hm$species_Eimeria <- gsub(pattern = "_", ". ",  hm$species_Eimeria)

hm$species_Eimeria <- gsub(pattern = "uninfected", "Uninfected",
                           hm$species_Eimeria)

hm$eimeriaSpecies <- gsub(pattern = "_", ". ",  hm$eimeriaSpecies)

hm$eimeriaSpecies <- gsub(pattern = "uninfected", "Uninfected",
                          hm$eimeriaSpecies)

# Then, define the color for each level of infection
color_mapping <- c("Uninfected controls" = "skyblue",
                   "E. ferrisi" = "forestgreen",
                   "E. falciformis" = "salmon")

factor_levels <- c("Uninfected controls",
                   "E. ferrisi",
                   "E. falciformis")

## repeat for field infections
# Then, define the color for each level of infection
color_mapping_f <- c("Uninfected" = "skyblue",
                   "E. ferrisi" = "forestgreen",
                   "E. falciformis" = "salmon")

factor_levels_f <- c("Uninfected",
                   "E. ferrisi",
                   "E. falciformis")

## Create factor levels for variables
Challenge$Parasite_primary <- 
    factor(Challenge$Parasite_primary, 
           levels = factor_levels)


Challenge$Parasite_challenge <- 
    factor(Challenge$Parasite_challenge, 
           levels = factor_levels)

hm$current_infection <- 
    factor(hm$current_infection, 
           levels = factor_levels)

hm$eimeriaSpecies <- 
    factor(hm$eimeriaSpecies, 
           levels = factor_levels_f)

hm$species_Eimeria <- 
    factor(hm$species_Eimeria, 
           levels = factor_levels_f)

# Change the species to italics in plot legends
labels = c("Uninfected controls",
             "*E. ferrisi*",
             "*E. falciformis*")

labels_f = c("Uninfected",
                    "*E. ferrisi*",
                    "*E. falciformis*")


### Create factors in mouse strains
# transform mouse strain into factor
hm$mouse_strain <- as.factor(hm$mouse_strain)

# Set NMRI as the reference level for mouse_strain
hm$mouse_strain <- relevel(hm$mouse_strain, ref = "NMRI")


# create lab
lab <- hm %>% 
    filter(origin == "Lab")

# update the immunization to truly reflect the uninfected mice
hm <- hm %>%
    mutate(immunization = case_when(
        Parasite_primary == "Uninfected controls" & Parasite_challenge == 
            "Uninfected controls" ~ "Uninfected controls",
        TRUE ~ coalesce(immunization, NA_character_)
    ))




## set the factor levels for the immunization variable
hm$immunization <- as.factor(hm$immunization)
hm$immunization <- relevel(hm$immunization, ref = "Uninfected controls")

