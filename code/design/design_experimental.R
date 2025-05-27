# Load your data (assuming lab_clean is your laboratory data)
library(dplyr)

# Basic infection counts from Table 2 data
infection_counts <- lab_clean %>%
    count(current_infection) %>%
    arrange(current_infection)

print("Basic infection counts:")
print(infection_counts)

# Get weight loss statistics by infection type
weight_loss_stats <- lab_clean %>%
    group_by(current_infection) %>%
    summarise(
        n = n(),
        mean_WL = mean(WL_max, na.rm = TRUE),
        min_WL = min(WL_max, na.rm = TRUE),
        max_WL = max(WL_max, na.rm = TRUE),
        sd_WL = sd(WL_max, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))

print("Weight loss statistics by infection:")
print(weight_loss_stats)

# Get weight loss by infection history (primary vs challenge)
# First, let's see what infection_history values you have
unique_infection_history <- lab_clean %>%
    distinct(infection_history, immunization) %>%
    arrange(infection_history)

print("Unique infection histories:")
print(unique_infection_history)

# Get weight loss by primary vs challenge infections
primary_challenge_stats <- lab_clean %>%
    filter(!is.na(immunization)) %>%
    group_by(current_infection, immunization) %>%
    summarise(
        n = n(),
        mean_WL = mean(WL_max, na.rm = TRUE),
        min_WL = min(WL_max, na.rm = TRUE),
        max_WL = max(WL_max, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))

print("Weight loss by infection type and immunization status:")
print(primary_challenge_stats)

# Get the specific numbers you need for your text
# Total mice infected
total_ferrisi <- sum(lab_clean$current_infection == "E. ferrisi", na.rm = TRUE)
total_falciformis <- sum(lab_clean$current_infection == "E. falciformis", na.rm = TRUE)
total_uninfected <- sum(lab_clean$current_infection == "Uninfected controls", na.rm = TRUE)

cat("\n=== NUMBERS FOR YOUR RESULTS TEXT ===\n")
cat("Total E. ferrisi infected:", total_ferrisi, "\n")
cat("Total E. falciformis infected:", total_falciformis, "\n") 
cat("Total uninfected controls:", total_uninfected, "\n")

# Get min/max for each group
ferrisi_stats <- lab_clean %>%
    filter(current_infection == "E. ferrisi") %>%
    summarise(
        mean = round(mean(WL_max, na.rm = TRUE), 1),
        min = round(min(WL_max, na.rm = TRUE), 1),
        max = round(max(WL_max, na.rm = TRUE), 1)
    )

falciformis_stats <- lab_clean %>%
    filter(current_infection == "E. falciformis") %>%
    summarise(
        mean = round(mean(WL_max, na.rm = TRUE), 1),
        min = round(min(WL_max, na.rm = TRUE), 1),
        max = round(max(WL_max, na.rm = TRUE), 1)
    )

uninfected_stats <- lab_clean %>%
    filter(current_infection == "Uninfected controls") %>%
    summarise(
        mean = round(mean(WL_max, na.rm = TRUE), 1),
        min = round(min(WL_max, na.rm = TRUE), 1),
        max = round(max(WL_max, na.rm = TRUE), 1)
    )

cat("\nE. ferrisi stats - Mean:", ferrisi_stats$mean, "%, Min:", ferrisi_stats$min, "%, Max:", ferrisi_stats$max, "%\n")
cat("E. falciformis stats - Mean:", falciformis_stats$mean, "%, Min:", falciformis_stats$min, "%, Max:", falciformis_stats$max, "%\n")
cat("Uninfected stats - Mean:", uninfected_stats$mean, "%, Min:", uninfected_stats$min, "%, Max:", uninfected_stats$max, "%\n")

# Check if you have separate primary and challenge data
if("Parasite_primary" %in% colnames(lab_clean)) {
    primary_stats <- lab_clean %>%
        filter(!is.na(Parasite_primary), Parasite_primary != "uninfected") %>%
        group_by(Parasite_primary) %>%
        summarise(
            n = n(),
            mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
            min_WL = round(min(WL_max, na.rm = TRUE), 1),
            max_WL = round(max(WL_max, na.rm = TRUE), 1),
            .groups = "drop"
        )
    
    cat("\nPRIMARY INFECTION STATS:\n")
    print(primary_stats)
}

if("Parasite_challenge" %in% colnames(lab_clean)) {
    challenge_stats <- lab_clean %>%
        filter(!is.na(Parasite_challenge), Parasite_challenge != "uninfected") %>%
        group_by(Parasite_challenge) %>%
        summarise(
            n = n(),
            mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
            min_WL = round(min(WL_max, na.rm = TRUE), 1),
            max_WL = round(max(WL_max, na.rm = TRUE), 1),
            .groups = "drop"
        )
    
    cat("\nCHALLENGE INFECTION STATS:\n")
    print(challenge_stats)
}

# Complete analysis to understand the experimental design
library(dplyr)

# 1. Basic overview
cat("=== BASIC OVERVIEW ===\n")
cat("Total mice:", nrow(lab_clean), "\n")
table(lab_clean$current_infection)

# 2. Understand the experimental structure
cat("\n=== EXPERIMENTAL STRUCTURE ===\n")
cat("Unique infection histories:\n")
table(lab_clean$infection_history)

cat("\nImmunization status:\n") 
table(lab_clean$immunization)

# 3. Cross-tabulation to see the full picture
cat("\n=== INFECTION HISTORY vs CURRENT INFECTION ===\n")
cross_tab <- table(lab_clean$infection_history, lab_clean$current_infection)
print(cross_tab)

# 4. Primary vs Challenge breakdown
cat("\n=== PRIMARY INFECTION BREAKDOWN ===\n")
if("Parasite_primary" %in% colnames(lab_clean)) {
    primary_breakdown <- lab_clean %>%
        count(Parasite_primary, current_infection) %>%
        arrange(Parasite_primary, current_infection)
    print(primary_breakdown)
}

cat("\n=== CHALLENGE INFECTION BREAKDOWN ===\n")
if("Parasite_challenge" %in% colnames(lab_clean)) {
    challenge_breakdown <- lab_clean %>%
        count(Parasite_challenge, current_infection) %>%
        arrange(Parasite_challenge, current_infection)
    print(challenge_breakdown)
}

# 5. Check for mice that had primary only vs primary + challenge
cat("\n=== PRIMARY vs PRIMARY+CHALLENGE ===\n")
primary_challenge_status <- lab_clean %>%
    mutate(
        had_primary = !is.na(Parasite_primary) & Parasite_primary != "uninfected",
        had_challenge = !is.na(Parasite_challenge) & Parasite_challenge != "uninfected",
        experiment_type = case_when(
            had_primary & had_challenge ~ "Primary + Challenge",
            had_primary & !had_challenge ~ "Primary only",
            !had_primary & had_challenge ~ "Challenge only", 
            TRUE ~ "Neither"
        )
    ) %>%
    count(experiment_type)

print(primary_challenge_status)

# 6. Weight loss by the ACTUAL experimental groups
cat("\n=== WEIGHT LOSS BY ACTUAL EXPERIMENTAL DESIGN ===\n")
weight_by_design <- lab_clean %>%
    group_by(infection_history, current_infection) %>%
    summarise(
        n = n(),
        mean_WL = round(mean(WL_max, na.rm = TRUE), 1),
        min_WL = round(min(WL_max, na.rm = TRUE), 1),
        max_WL = round(max(WL_max, na.rm = TRUE), 1),
        .groups = "drop"
    ) %>%
    arrange(current_infection, infection_history)

print(weight_by_design)


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
###################################################











        










# Combine the plots
#(ooc_primary | ooc_challenge) / # oocysts
panel_figure <- 
    (Rwp | Rwc) /
    (strains_weight_challenge ) / 
    (eimeria_weight) +
 #   plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure <- panel_figure + 
    plot_annotation(title = 'Fig. 1', 
                    theme = theme(plot.title = element_text(size = 20,
                                                            hjust = 0)))

# Control sizes of each plot within the panel
# This is a generic example. You'll need to adjust the widths,
#heights, and layout design based on your specific needs.
#panel_figure <- panel_figure + 
 #   plot_layout(heights = c(1, 1, 1), 
  #              widths = c(1, 1, 1)) # Adjust according to your layout needs

# Display the panel figure
print(panel_figure)

# Save the panel figure
ggsave(paste0(panels_fi, "/experimental_design_simple.jpeg"), 
       panel_figure, width = 10, height = 12, dpi = 300)

############################
#Hybrids


# Combine the plots
panel_hybr <- 
    (h_w | map_hybrids) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_hybr <- panel_hybr + 
    plot_annotation(title = 'Fig. 7', 
                    theme = theme(plot.title = 
                                      element_text(size = 13, hjust = 0)))

# Display the panel figure
print(panel_hybr)

# Save the panel figure
ggsave(paste0(panels_fi, "/hyb_plot_design.jpeg"), 
       panel_hybr, width = 6, height = 3, dpi = 300)



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



