# Filter the data for relevant columns and origin
wild_and_lab <- hm %>%
    dplyr::ungroup() %>%
    dplyr::select(Mouse_ID, origin, infection_status, species_Eimeria, all_of(Genes_v))  # Include infection variables
# Include genes of interest

# Perform PCA on the gene expression data (excluding 'origin' column for PCA)
pca_input <- wild_and_lab %>%
    ungroup() %>%
    dplyr::select(all_of(Genes_v))

# Run PCA
res_pca <- PCA(pca_input, graph = FALSE)

# Extract PCA coordinates
pca_coordinates <- as.data.frame(res_pca$ind$coord) %>%
    dplyr::mutate(origin = wild_and_lab$origin, 
                  infection_status = wild_and_lab$infection_status,
                  species_Eimeria = wild_and_lab$species_Eimeria)  # Add origin for grouping

# Modify PCA plot
pca_plot <- ggplot(pca_coordinates, aes(x = Dim.1, y = Dim.2, color = origin, fill = origin)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.3, linetype = "dashed") +
    labs(
        x = paste0("PC1 (", round(res_pca$eig[1, 2], 1), "% variance explained)"),
        y = paste0("PC2 (", round(res_pca$eig[2, 2], 1), "% variance explained)"),
        color = "Sample Type",  # Rename "Origin"
        fill = "Sample Type"    # Rename "Origin"
    ) +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Use a beautiful viridis palette
    scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +   # Use the same palette for fill
    theme_minimal(base_size = 16) +
    theme(
        legend.position = "top",
        plot.title = element_blank(),  # Remove the title
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
    )

# Display the modified plot
print(pca_plot)

# Save the modified plot
ggsave(
    filename = paste0(an_fi, "/PCA_Wild_Lab_Gene_Expression_Modified.jpeg"),
    plot = pca_plot,
    width = 8,
    height = 6,
    dpi = 300
)




###############################################################################
# creating variables showing origin and infection status
# Create a new variable combining 'origin' and 'infection_status'
pca_coordinates <- pca_coordinates %>%
    dplyr::mutate(
        origin_infection = case_when(
            infection_status == TRUE ~ paste0(origin, "_infected"),
            infection_status == FALSE ~ paste0(origin, "_uninfected")
        ), 
        origin_species = paste0(origin, species_Eimeria)
    )

# Check the updated dataframe
head(pca_coordinates)


# Modify PCA plot with shape based on infection_status
ggplot(pca_coordinates %>%
           drop_na(infection_status), aes(x = Dim.1, y = Dim.2, color = origin_infection, fill = origin_infection)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.3, linetype = "dashed") +
    labs(
        x = paste0("PC1 (", round(res_pca$eig[1, 2], 1), "% variance explained)"),
        y = paste0("PC2 (", round(res_pca$eig[2, 2], 1), "% variance explained)"),
        color = "Sample Type",  # Rename "Origin"
        shape = "Infection Status",  # Add legend for shape
        fill = "Sample Type"
    ) +
    #scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Use a beautiful viridis palette
    #scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +   # Use the same palette for fill
    scale_shape_manual(values = c(16, 17)) +  # Map shapes (e.g., circle = 16, triangle = 17)
    theme_minimal(base_size = 16) +
    theme(
        legend.position = "top",
        plot.title = element_blank(),  # Remove the title
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
    )



# Modify PCA plot with shape based on infection_status
ggplot(pca_coordinates %>%
           drop_na(species_Eimeria), aes(x = Dim.1, y = Dim.2, color = origin_species, fill = origin_species)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.3, linetype = "dashed") +
    labs(
        x = paste0("PC1 (", round(res_pca$eig[1, 2], 1), "% variance explained)"),
        y = paste0("PC2 (", round(res_pca$eig[2, 2], 1), "% variance explained)"),
        color = "Sample Type",  # Rename "Origin"
        shape = "Infection Status",  # Add legend for shape
        fill = "Sample Type"
    ) +
    #scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Use a beautiful viridis palette
    #scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +   # Use the same palette for fill
    scale_shape_manual(values = c(16, 17)) +  # Map shapes (e.g., circle = 16, triangle = 17)
    theme_minimal(base_size = 16) +
    theme(
        legend.position = "top",
        plot.title = element_blank(),  # Remove the title
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
    )




