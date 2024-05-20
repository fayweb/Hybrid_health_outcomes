# 7.4: Here we test the hybrid impact on predicted weight loss and further test 
# the effect of sex and infection status with Eimeria 
# I have further tested to see if the presence of various parasites independent 
# of eimeria have an impact on the predicted weight loss in combination with 
# hybridicity 


# In case there are issues with the package ParasiteLoad, you may have to update
# the package optimx to the latter version
#install.packages("optimx", version = "2021-10.12")

## install the pacakage of Alice Balard
#devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)
#force = TRUE)


# Testing differences between female and male hybrids of M.m. musculus and 
# m.m.domesticus in predicted weight loss
Field$Sex <- as.factor(Field$Sex)


##All
fitWL_Sex <- parasiteLoad::analyse(data = Field,
                                   response = "predicted_WL",
                                   model = "normal",
                                   group = "Sex")
     

plot_WL_Sex<- bananaPlot(mod = fitWL_Sex$H3,
                         data = Field,
                         response = "predicted_WL",
                         group = "Sex",
                         cols = c("white", "white")) +
    scale_fill_manual(values = c("orange", "forestgreen"),
                      name = "Sex") +
    scale_color_manual(values = c("orange", "forestgreen"),
                       name = "Sex") +
    theme_bw() +
    labs(y = "Predicted weight loss") 

plot_WL_Sex


### create the hybrid bar
# Adjust the gradient bar plot to include axis labels and remove space
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)), 
                        aes(x=hi, y=1, fill = hi)) +
    geom_tile() +
    scale_x_continuous(breaks=seq(0, 1, by=0.25), 
                       labels=c("0", "0.25", "0.5", "0.75", "1")) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_void() +
    theme(legend.position = 'none', 
          plot.margin = unit(c(-1, 0, 0, 0), "npc"), 
          # This removes space around the plot
          axis.text.x = element_text(color = "black", 
                                     angle = 0, vjust = 0.5, hjust=0.5)) 
# Adjust text vjust for positioning

HIgradientBar

plot_WL_Sex <- plot_WL_Sex +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))


plot_WL_Sex

#
ggsave(plot = plot_WL_Sex, filename = paste0(an_fi, "/hybrid_sex.jpeg"), 
       width = 6, 
       height = 6, dpi = 1000)

# Ensure the bottom plot (HIgradientBar) is ready (assuming HIgradientBar 
#is already defined)
HIgradientBar <- HIgradientBar + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Use patchwork to combine the plots without any space between them
combined_plot <- plot_WL_Sex / HIgradientBar + 
    plot_layout(heights = c(1, 0.1)) # Adjust the relative heights as needed

# Print the combined plot
combined_plot


ggsave(plot = combined_plot, 
       filename = paste0(an_fi, "/hybrid_sex.jpeg"), width = 4, height = 4, 
       dpi = 1000)


####################### Mapping ######################################
leaflet(data = Field) %>%
    addTiles() %>%
    addMarkers(
        lng = ~Longitude, lat = ~Latitude, popup = ~as.character(Mouse_ID))

colorPalette <- colorRampPalette(c("blue", "red"))


colors <- colorPalette(100)[as.numeric(cut(Field$HI, breaks = 100))]

leaflet_map <-
    leaflet(data = Field) %>%
    addTiles() %>%
    addCircleMarkers(lng = ~Longitude, lat = ~Latitude, color = ~colors, 
                     radius = 5, fillOpacity = 0.8, stroke = FALSE, 
                     popup = ~as.character(HI))


#Read the Leaflet map image
#leaflet_image <- magick::image_read(paste0(an_fi, "/Hybrid_map.jpeg"))

# Convert to a raster for grid plotting
#leaflet_raster <- rasterGrob(leaflet_image, interpolate = TRUE)

#leaflet_raster <- ggplot() +
 #   background_image(leaflet_image) + coord_fixed()
###################################



##############################################################################
### Testing diffences between infected and uninfected hybrid mice

# Melting Curve analysis
Field_mc <- Field %>%
    drop_na(MC.Eimeria)


Field_mc <- Field_mc %>%
    mutate(MC.Eimeria = if_else(MC.Eimeria == "TRUE", "infected", "uninfected"))

Field_mc$MC.Eimeria <- as.factor(Field_mc$MC.Eimeria)
fitWL_mc <- parasiteLoad::analyse(data = Field_mc,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "MC.Eimeria")


plot_WL_mc <- 
    bananaPlot(mod = fitWL_mc$H3,
               data = Field_mc,
               response = "predicted_WL",
               group = "MC.Eimeria",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("purple", "cornflowerblue"), 
                      name = "Eimeria spp.") +
    scale_color_manual(values = c("purple", "cornflowerblue"),
                       name = "Eimeria spp.") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted weight loss")

plot_WL_mc

# Create the combined plot with the gradient bar as the "axis"
plot_WL_mc_combined <-  (plot_WL_mc / HIgradientBar) + 
    plot_layout(heights = c(1, 0.1)) # Adjust the relative heights as needed


# Display the combined plot
#plot_WL_mc_combined


ggsave(plot = plot_WL_mc_combined, 
       filename = paste0(an_fi, "/hybrid_meltingcurve.jpeg"),  
       width = 6, height = 6, dpi = 1000)

Field <- Field %>%
    mutate(HE = 2*HI*(1-HI)) 

ggplot(Field, aes(x = HE, predicted_WL, color = Sex)) +
    geom_smooth(method = lm, se = TRUE) 

ggplot(Field, aes(x = HE, predicted_WL, color = infection_status)) +
    geom_smooth(method = lm, se = TRUE) 

######################
# Combine the plots without any space between them and add labels 
#only to the first row
map_plot <- (combined_plot | plot_WL_mc_combined) + 
    plot_layout(tag_level = 'new') +  # 
    plot_annotation(tag_levels = list(c('A', '', 'B', '')))
# +  # Add labels (A, B)

# Add a figure title
map_plot <- map_plot + 
    plot_annotation(title = 'Fig. 8', 
                    theme = theme(plot.title = element_text(size = 13, 
                                                            hjust = 0),
                                  legend.position = "none"))

# Display the panel figure
print(map_plot)

# Save the plot
# Note: Corrected the paste0 function usage for filename specification
ggsave(plot = map_plot, 
       filename = paste0(panels_fi, "/banana_map_immune_signature.jpeg"), 
       width = 10, 
       height = 5, 
       dpi = 1000)

#############################
########################
###################
################ combining amplicon data
### Testing diffences between infected and uninfected hybrid mice

# Melting Curve analysis
Field_mc <- Field %>%
    drop_na(infection_status)


Field_mc <- Field_mc %>%
    mutate(infection_status = if_else(infection_status == "TRUE", "infected", "uninfected"))

Field_mc$infection_status <- as.factor(Field_mc$infection_status)

fitWL_mc <- parasiteLoad::analyse(data = Field_mc,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "infection_status")


plot_WL_mc <- 
    bananaPlot(mod = fitWL_mc$H3,
               data = Field_mc,
               response = "predicted_WL",
               group = "infection_status",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("purple", "cornflowerblue"), 
                      name = "Eimeria spp.") +
    scale_color_manual(values = c("purple", "cornflowerblue"),
                       name = "Eimeria spp.") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted weight loss")

plot_WL_mc

# Create the combined plot with the gradient bar as the "axis"
plot_WL_mc_combined <-  (plot_WL_mc / HIgradientBar) + 
    plot_layout(heights = c(1, 0.1)) # Adjust the relative heights as needed


# Display the combined plot
plot_WL_mc_combined


ggsave(plot = plot_WL_mc_combined, 
       filename = paste0(an_fi, "/hybrid_meltingcurve.jpeg"),  
       width = 6, height = 6, dpi = 1000)


######################
# Combine the plots without any space between them and add labels 
#only to the first row
map_plot <- (combined_plot | plot_WL_mc_combined) + 
    plot_layout(tag_level = 'new') +  # 
    plot_annotation(tag_levels = list(c('A', '', 'B', '')))
# +  # Add labels (A, B)

# Add a figure title
map_plot <- map_plot + 
    plot_annotation(title = 'Fig. 8', 
                    theme = theme(plot.title = element_text(size = 13, 
                                                            hjust = 0),
                                  legend.position = "none"))

# Display the panel figure
print(map_plot)

# Save the plot
# Note: Corrected the paste0 function usage for filename specification
ggsave(plot = map_plot, 
       filename = paste0(panels_fi, "/banana_map_immune_signature.jpeg"), 
       width = 10, 
       height = 5, 
       dpi = 1000)






##################################################       
############ Testing according to delta ct (not recommended, just validating)
Field_ct <-  Field %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    dplyr::mutate(delta_infection = 
                      case_when(delta_ct_cewe_MminusE < -5 ~ "uninfected",
                                delta_ct_cewe_MminusE > -5 ~ "infected"))

Field_ct$delta_infection <- as.factor(Field_ct$delta_infection)

##All
fitWL_ct<- parasiteLoad::analyse(data = Field_ct,
                                 response = "predicted_WL",
                                 model = "normal",
                                 group = "delta_infection")


plot_WL_ct <- 
    bananaPlot(mod = fitWL_ct$H3,
               data = Field_ct,
               response = "predicted_WL",
               group = "delta_infection",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("indianred3","steelblue1"),
                      name = "qPCR values - Eimeria detection in Caecum") +
    scale_color_manual(values = c("indianred3", "steelblue1"),
                       name = "qPCR values - Eimeria detection in Caecum") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_ct

combined_plot_ct <- plot_WL_ct / 
    HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))


ggsave(plot = combined_plot_ct, 
       filename = paste0(an_fi, "/hybrid_infected_ct.jpeg"),  
       width = 6, height = 6, dpi = 1000)


##################### Eimeria oocysts
Field_ooc_eim <-  Field %>%
    drop_na(OPG) %>%
    dplyr::mutate(oocysts_present = 
                      case_when(OPG > 0 ~ "TRUE",
                                OPG == 0 ~ "FALSE"))

Field_ooc_eim$oocysts_present <- as.factor(Field_ooc_eim$oocysts_present)

##All
fitWL_ooc<- parasiteLoad::analyse(data = Field_ooc_eim,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "oocysts_present")

plot_WL_ooc <- 
    bananaPlot(mod = fitWL_ooc$H3,
               data = Field_ooc_eim,
               response = "predicted_WL",
               group = "oocysts_present",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1","indianred3"),
                      name = "Eimeria oocysts present") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Eimeria oocysts present") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_ooc


# Create the combined plot with the gradient bar as the "axis"
combined_plot_ooc <- plot_WL_ooc / 
    HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
combined_plot_ooc


ggsave(plot = combined_plot_ooc, 
       filename = paste0(an_fi, "/hybrid_eimeria_ooc.jpeg"),
       width = 6, height = 6, 
       dpi = 1000)

###################################

###################### I will here try to test it with other parsite infection
## I use different ways to code the presence of absence with the parasite
# note that this doesn't use a quantification measure
####################################################
#####  crpyto ct
Field_crypto <- Field %>%
    drop_na(ILWE_Crypto_Ct) %>%
    mutate(crypto_infected_ct = 
               case_when(
                   ILWE_Crypto_Ct > 0 ~ "TRUE",
                   ILWE_Crypto_Ct == 0 ~ "FALSE",
                   ILWE_Crypto_Ct == NA ~ NA
               ))

Field_crypto <- Field_crypto %>%
    drop_na(crypto_infected_ct)

Field_crypto$crypto_infected_ct <- as.factor(Field_crypto$crypto_infected_ct)



fitWL_crypto_ct <- parasiteLoad::analyse(data = Field_crypto,
                                         response = "predicted_WL",
                                         model = "normal",
                                         group = "crypto_infected_ct")


# plot it
plot_WL_crypto_ct <- 
    bananaPlot(mod = fitWL_crypto_ct$H3,
               data = Field_crypto,
               response = "predicted_WL",
               group = "crypto_infected_ct",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Cryptosporidum detected in qPCR") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Cryptosporidum detected in qPCR") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_crypto_ct

# Create the combined plot with the gradient bar as the "axis"
plot_WL_crypto_Ct <- plot_WL_crypto_ct / 
    HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
plot_WL_crypto_ct


ggsave(plot = plot_WL_crypto_Ct, 
       filename = paste0(an_fi, "/hybrid_crypto_qpcr.jpeg"),  width = 6, 
       height = 6, dpi = 1000)


################################################
##############
######Worms?
##################First aspiculuris
Field %>%
    pivot_longer(cols = c("Trichuris_muris", "Mastophorus_muris", 
                          "Catenotaenia_pusilla", "Heligmosomoides_polygurus",
                          "Heterakis_sp", "Aspiculuris_sp", "Syphacia_sp",
                          "Taenia_sp", "Hymenolepis_sp"), names_to = "Worms", 
                 values_to = "Counts_worms") %>%
    ggplot(aes(x = Worms, y = Counts_worms)) +
    geom_violin()




#####  Aspiculuris_sp
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Aspiculuris_sp > 0 ~ "TRUE",
                   Aspiculuris_sp == 0 ~ "FALSE",
                   Aspiculuris_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_worms <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Aspiculuris") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Aspiculuris") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_worms

# Create the combined plot with the gradient bar as the "axis"
plot_WL_worms <- plot_WL_worms / HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
plot_WL_worms


ggsave(plot = plot_WL_worms, 
       filename = paste0(an_fi, "/hybrid_aspiculuris.jpeg"),  
       width = 6, height = 6, 
       dpi = 1000)

##############################################################
#####  Syphacia
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Syphacia_sp > 0 ~ "TRUE",
                   Syphacia_sp == 0 ~ "FALSE",
                   Syphacia_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_syphacia <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Syphacia") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Syphacia") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_syphacia

# Create the combined plot with the gradient bar as the "axis"
plot_WL_syphacia <-  plot_WL_syphacia / HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
plot_WL_syphacia


ggsave(plot = plot_WL_syphacia, 
       filename = paste0(an_fi, "/hybrid_syphacia.jpeg"),
       width = 6, height = 6, 
       dpi = 1000)

##############  Trichuris
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Trichuris_muris > 0 ~ "TRUE",
                   Trichuris_muris == 0 ~ "FALSE",
                   Trichuris_muris == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")

# plot it
plot_WL_trichuris <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Trichuris") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Trichuris") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_trichuris

# Create the combined plot with the gradient bar as the "axis"
plot_WL_trichuris <- plot_WL_trichuris / HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
plot_WL_trichuris


ggsave(plot = plot_WL_trichuris, 
       filename = paste0(an_fi, "/hybrid_trichuris.jpeg"),
       width = 6, height = 6, dpi = 1000)


##############  Heterakis
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Heterakis_sp > 0 ~ "TRUE",
                   Heterakis_sp == 0 ~ "FALSE",
                   Heterakis_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_heterakis <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Heterakis") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Heterakis") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_heterakis

# Create the combined plot with the gradient bar as the "axis"
plot_WL_heterakis <- plot_WL_heterakis / HIgradientBar + 
    plot_layout(heights = c(1.3, 1/8))

# Display the combined plot
plot_WL_heterakis


ggsave(plot = plot_WL_heterakis, 
       filename = paste0(an_fi, "/hybrid_heterakis.jpeg"),
       width = 6, height = 6, dpi = 1000)



infection_plots <- 
    cowplot::plot_grid(plot_WL_mc_combined, combined_plot, combined_plot_ooc, 
              plot_WL_crypto_Ct, plot_WL_worms, plot_WL_syphacia, 
              plot_WL_trichuris,
              plot_WL_heterakis,
              nrow = 4, align = c("hv"), labels = LETTERS[1:8])

ggsave(infection_plots, filename = paste0(panels_fi, "/infection_panel.jpeg"),
       width = 12, height = 24, dpi = 1000)

