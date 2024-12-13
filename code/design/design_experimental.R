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



