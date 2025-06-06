# =============================================================================
# SCRIPT PURPOSE: In-depth exploration of experimental design for lab infections
# =============================================================================
# This script performs:
# 1. Summary statistics on infection groups
# 2. Tabulations of infection history, immunization, and challenge status
# 3. Weight loss statistics across groups
# 4. Generation of summary figures (PCA, weight loss, HI maps)
# 5. Model fitting for weight loss predictors
# =============================================================================
# === Load cleaned datasets ===
Challenge <- read.csv("data/analysis/final/Challenge_ready_for_analysis.csv")
hm <- read.csv("data/analysis/final/hm_ready_for_analysis.csv")


lab <- hm %>%
    dplyr::filter(origin == "Lab")
lab_clean <- lab %>% ungroup()

cat("=== BASIC SAMPLE SIZES ===\n")
cat("Total mice:", nrow(lab_clean), "\n\n")

# Final infection status counts
infection_counts <- lab_clean %>%
    count(current_infection) %>%
    arrange(current_infection)

cat("Final infection status:\n")
print(infection_counts)

# =============================================================================
# 2. EXPERIMENTAL STRUCTURE OVERVIEW
# =============================================================================

cat("\n=== EXPERIMENTAL STRUCTURE ===\n")

# Infection histories
cat("Infection histories (primary_challenge):\n")
infection_history_table <- table(lab_clean$infection_history)
print(infection_history_table)

# Immunization categories
cat("\nImmunization categories:\n")
immunization_table <- table(lab_clean$immunization)
print(immunization_table)

# Cross-tabulation: infection history vs final status
cat("\nInfection history vs Final status:\n")
cross_tab <- table(lab_clean$infection_history, lab_clean$current_infection)
print(cross_tab)

# =============================================================================
# 3. PRIMARY vs CHALLENGE BREAKDOWN
# =============================================================================

cat("\n=== PRIMARY vs CHALLENGE BREAKDOWN ===\n")

# Primary infection breakdown
if("Parasite_primary" %in% colnames(lab_clean)) {
    cat("Primary infections by final status:\n")
    primary_breakdown <- lab_clean %>%
        count(Parasite_primary, current_infection) %>%
        arrange(Parasite_primary, current_infection)
    print(primary_breakdown)
}

# Challenge infection breakdown  
if("Parasite_challenge" %in% colnames(lab_clean)) {
    cat("\nChallenge infections by final status:\n")
    challenge_breakdown <- lab_clean %>%
        count(Parasite_challenge, current_infection) %>%
        arrange(Parasite_challenge, current_infection)
    print(challenge_breakdown)
}

# Experimental design verification
cat("\nExperimental design verification:\n")
primary_challenge_status <- lab_clean %>%
    mutate(
        had_primary = !is.na(Parasite_primary) & Parasite_primary != "uninfected",
        had_challenge = !is.na(Parasite_challenge) & Parasite_challenge != "uninfected",
        experiment_type = case_when(
            had_primary & had_challenge ~ "Primary + Challenge",
            had_primary & !had_challenge ~ "Primary only",
            !had_primary & had_challenge ~ "Challenge only", 
            TRUE ~ "Control/Neither"
        )
    ) %>%
    count(experiment_type)

print(primary_challenge_status)

# =============================================================================
# 4. WEIGHT LOSS STATISTICS
# =============================================================================

cat("\n=== WEIGHT LOSS STATISTICS ===\n")

# By final infection status
cat("Weight loss by FINAL infection status:\n")
final_weight_stats <- lab_clean %>%
    group_by(current_infection) %>%
    summarise(
        n = n(),
        mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL = round(min(WL_max, na.rm = TRUE), 1),
        max_WL = round(max(WL_max, na.rm = TRUE), 1),
        sd_WL = round(sd(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    )
print(final_weight_stats)

# By immunization status
cat("\nWeight loss by immunization category:\n")
immunization_weight_stats <- lab_clean %>%
    group_by(current_infection, immunization) %>%
    summarise(
        n = n(),
        mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL = round(min(WL_max, na.rm = TRUE), 1),
        max_WL = round(max(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    ) %>%
    arrange(current_infection, immunization)
print(immunization_weight_stats)

# By detailed infection history
cat("\nWeight loss by detailed infection history:\n")
detailed_weight_stats <- lab_clean %>%
    group_by(infection_history, current_infection) %>%
    summarise(
        n = n(),
        mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL = round(min(WL_max, na.rm = TRUE), 1),
        max_WL = round(max(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    ) %>%
    arrange(current_infection, infection_history)
print(detailed_weight_stats)

# =============================================================================
# 5. SUMMARY FOR RESULTS WRITING
# =============================================================================

cat("\n=== SUMMARY FOR RESULTS TEXT ===\n")
cat("Total mice: 136\n")
cat("E. ferrisi final status:", sum(lab_clean$current_infection == "E. ferrisi"), "\n")
cat("E. falciformis final status:", sum(lab_clean$current_infection == "E. falciformis"), "\n")
cat("Uninfected controls:", sum(lab_clean$current_infection == "Uninfected controls"), "\n")

# Key weight loss numbers
ferrisi_summary <- lab_clean %>% 
    filter(current_infection == "E. ferrisi") %>%
    summarise(mean = round(mean(WL_max, na.rm = TRUE), 1),
              min = round(min(WL_max, na.rm = TRUE), 1),
              max = round(max(WL_max, na.rm = TRUE), 1))

falciformis_summary <- lab_clean %>% 
    filter(current_infection == "E. falciformis") %>%
    summarise(mean = round(mean(WL_max, na.rm = TRUE), 1),
              min = round(min(WL_max, na.rm = TRUE), 1),
              max = round(max(WL_max, na.rm = TRUE), 1))

controls_summary <- lab_clean %>% 
    filter(current_infection == "Uninfected controls") %>%
    summarise(mean = round(mean(WL_max, na.rm = TRUE), 1),
              min = round(min(WL_max, na.rm = TRUE), 1),
              max = round(max(WL_max, na.rm = TRUE), 1))

cat("\nWeight loss summaries:\n")
cat("E. ferrisi: mean =", ferrisi_summary$mean, "%, range =", ferrisi_summary$min, "-", ferrisi_summary$max, "%\n")
cat("E. falciformis: mean =", falciformis_summary$mean, "%, range =", falciformis_summary$min, "-", falciformis_summary$max, "%\n")
cat("Controls: mean =", controls_summary$mean, "%, range =", controls_summary$min, "-", controls_summary$max, "%\n")


# Check for mice that died/were sacrificed during primary infections
cat("=== UNDERSTANDING SACRIFICED MICE ===\n")

# Check if there's a 'death' or 'sacrifice' column
if("death" %in% colnames(lab_clean)) {
    cat("Death/sacrifice information:\n")
    table(lab_clean$death, useNA = "always")
}

# Check experiment column for any patterns
if("experiment" %in% colnames(lab_clean)) {
    cat("Mice by experiment:\n")
    table(lab_clean$experiment, useNA = "always")
}

# Look at mice that might have had primary only
cat("\nDetailed breakdown of ALL 136 mice:\n")
detailed_breakdown <- lab_clean %>%
    count(Parasite_primary, Parasite_challenge, current_infection) %>%
    arrange(Parasite_primary, Parasite_challenge, current_infection)

print(detailed_breakdown)

# Check for any mice with missing challenge data
cat("\nMice with missing challenge infection data:\n")
missing_challenge <- lab_clean %>%
    filter(is.na(Parasite_challenge)) %>%
    count(Parasite_primary, current_infection)

if(nrow(missing_challenge) > 0) {
    print(missing_challenge)
} else {
    cat("No mice with missing challenge data found.\n")
}

# Total verification
cat("\nTotal mice verification:\n")
cat("Total mice in dataset:", nrow(lab_clean), "\n")
cat("Mice with primary infection data:", sum(!is.na(lab_clean$Parasite_primary)), "\n")
cat("Mice with challenge infection data:", sum(!is.na(lab_clean$Parasite_challenge)), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")


# Select laboratory data 
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

Field <- hm %>%
    dplyr::filter(origin == "Field")

# Change the species to italics in plot legends
labels = c("Uninfected controls",
              "*E. ferrisi*",
              "*E. falciformis*")
                 

lab$Parasite_primary <- 
    factor(lab$Parasite_primary, 
           levels = factor_levels)

lab$Parasite_challenge <- 
    factor(lab$Parasite_challenge, 
           levels = factor_levels)

# check the distributions of weight loss for each mouse strain
# Define colors
colors <- c("TRUE" = "firebrick3", "FALSE" = "steelblue")




# Creating a density plot for the Hybrid Index (HI)
ggplot(Field, aes(HI)) + 
    geom_density(fill = "steelblue",alpha = 0.7) +
    geom_vline(aes(xintercept = mean(HI, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    labs(#title = "Distribution of Hybrid Index (HI) Among Wild Mice",
        x = "Hybrid Index (HI)",
        y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) -> h_w

h_w


ggsave(filename =  paste0(d_fi,"/densityplot_HI.jpeg"),
       plot = h_w, width = 8, height = 6)
# The red dashed line represents the mean HI value, providing a reference for 
#the central tendency of hybridization.


# Base world map
world_map <- map_data("world")

lon_range <- range(Field$Longitude) + c(-0.1, 0.1)  # Expanding the range a 
#bit for padding
lat_range <- range(Field$Latitude) + c(-0.1, 0.1)   # Expanding the range a 
#bit for padding

# Plot with zoom
map_hybrids <-
    ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
                 fill = "gray90", color = "white") +
    geom_point(data = Field, aes(x = Longitude, y = Latitude, color = HI), 
               size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(#title = "Mouse HI Values Across Locations", 
        x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_fixed(1.3) +  # Keep map aspect ratio
    xlim(lon_range) +   # Set x-axis limits to zoom in
    ylim(lat_range)     # Set y-axis limits to zoom in

map_hybrids 

ggsave(filename = paste0(d_fi,"/map_HI.jpeg"),
       plot = map_hybrids, width = 8, height = 6)

###################################################
####################################################
##################################################





# prim 
#prim <- 
#count(lab$Parasite_primary)


#########################################################################
Challenge <- Challenge %>%
    group_by(Mouse_ID, infection) %>%
    mutate(weight_loss = 100 - min(relative_weight))

## statistics primary infection
Challenge %>%
    filter(infection == "primary") %>%
    group_by(Parasite_primary) %>%
    summarise(
        MeanWeightLoss = mean(weight_loss),
        MinWeightLoss = min(weight_loss),
        MaxWeightLoss = max(weight_loss),
        N = n()
    )

## statistics challenge infection
Challenge %>%
    filter(infection == "challenge") %>%
    group_by(Parasite_challenge) %>%
    summarise(
        MeanWeightLoss = mean(weight_loss),
        MinWeightLoss = min(weight_loss),
        MaxWeightLoss = max(weight_loss),
        N = n()
    )

## statistics primary infection - getting the N
lab %>%
    filter(infection == "primary") %>%
    group_by(Parasite_primary) %>%
    summarise(
        MeanWeightLoss = mean(WL_max),
        MinWeightLoss = min(WL_max),
        MaxWeightLoss = max(WL_max),
        N = n()
    )

# - getting the N
lab %>%
    filter(infection == "challenge")%>%
    group_by(Parasite_challenge) %>%
    summarise(
        MeanWeightLoss = mean(WL_max),
        MinWeightLoss = min(WL_max),
        MaxWeightLoss = max(WL_max),
        N = n()
    )

# mean dpi at peak weight loss falciformis
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "primary", Parasite_primary == "E_falciformis",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss falciformis - challenge
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "challenge", Parasite_challenge == "E_falciformis",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss primary ferrisi
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "primary", Parasite_primary == "E_ferrisi",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss E_ferrisi - challenge
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "challenge", Parasite_challenge == "E_ferrisi",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

######################models
###############

model <- lm(WL_max ~ mouse_strain, data = lab)
summary(model)

# create a new combined variable 
Challenge_p <- Challenge %>%
    filter(infection == "primary") %>%
    mutate(Parasite = Parasite_primary)

Challenge_c <- Challenge %>%
    filter(infection == "challenge") %>%
    mutate(Parasite = Parasite_challenge)

Challenge <- rbind(Challenge_p, Challenge_c)
rm(Challenge_c, Challenge_p)

model2 <- lm(weight_loss ~ infection * Parasite, data = Challenge)
summary(model2)




rm(chale, challenge, Eim_strains, eimeria_weight_chal, eimeria_weight_challenge,
   eimeria_weight_prim, h_w, m_s, map_hybrids, model, model2, mouse_WL, 
   ooc_challenge, ooc_primary, panel_figure, panel_hybr, parasite_WL, primary, 
   Rwc, Rwp, s, strains_weight, strains_weight_challenge, world_map)



