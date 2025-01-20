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


#################################################   Plot 1A
#########################################
# Let's add both challenge and primary together 
lab  %>% 
    ggplot(aes(x = WL_max, y = current_infection, fill = current_infection)) + 
    geom_density_ridges(jittered_points = TRUE, 
                        position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, 
                        point_size = 2, point_alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7,
                 position = position_nudge(x = 0.1)) +
    theme_minimal(base_size = 16) +
    scale_fill_manual(values = color_mapping) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 20, face = "bold"),  # X-axis title size
          axis.title.y = element_text(size = 20, face = "bold"),  # Y-axis title size
          axis.text.x = element_text(size = 20),  # X-axis tick size
          axis.text.y = element_text(size = 20))+ # Y-axis tick size) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain")  -> eimeria_weight

eimeria_weight

eimeria_weight <- italics_y(eimeria_weight, labels)
eimeria_weight

ggsave(filename = paste0(d_fi,"/eimeria_weight_combined.pdf"),
       plot = eimeria_weight, 
       width = 8, height = 6, dpi = 300)

###############################################      Plot 1B
lab <- hm %>%
    filter(origin == "Lab")

# select the laboratory genes
genes <- lab[ ,colnames(lab) %in% Genes_v]

# PCA
##pca on the complete data set
res.pca <- PCA(genes)

############### biplot
biplot <- fviz_pca_biplot(
    res.pca, 
    col.ind = lab$current_infection,
    pointsize = 2,
    addEllipses = TRUE, 
    ellipse.level = 0.8, # Set ellipse confidence level to 95%
    alpha.ind = 0.8,
    alpha.var = 1, 
    label = "var",
    col.var = "black", 
    repel = TRUE,
    legend.title = "Parasite strain",
    title = ""
) +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    scale_shape_manual(values = c(15, 16, 17), labels = labels) +
    labs(color = "Parasite strain", shape = "Parasite strain") +
    theme_minimal(base_size = 16) +
    theme(legend.text = element_markdown(),
          axis.title.x = element_text(size = 20, face = "bold"),  # X-axis title size
          axis.title.y = element_text(size = 20, face = "bold"),  # Y-axis title size
          axis.text.x = element_text(size = 20),  # X-axis tick size
          axis.text.y = element_text(size = 20)) # Enable Markdown for italics

biplot

ggsave(filename = paste0(an_fi, "/biplot.pdf"), plot = biplot, 
       width = 8, height = 6, dpi = 300)


# coefficient Plot 1 
## make a coefficient plot with just the PC1 and PC2
# only pc1 + pc2
# Convert mouse_id to a data frame
mouse <- data.frame(Mouse_ID = lab[,1])

# Add the new column pc1 to the mouse_id data frame
mouse$PC1 <- res.pca$ind$coord[, 1]
mouse$PC2 <- res.pca$ind$coord[, 2]  # indexing the second column
mouse$PC3 <-  res.pca$ind$coord[, 3]
mouse$PC4 <-  res.pca$ind$coord[, 4]
mouse$PC5 <-  res.pca$ind$coord[, 5]

lab <- lab %>% 
    left_join(mouse, by = "Mouse_ID")





#################### Interaction effects model
model_6 <- lm(WL_max ~ PC1 * current_infection + PC2 *current_infection, 
              data = lab)

summary(model_6)


# Set consistent axis limits
x_limits <- c(-8, 5)  # Example limits for the x-axis (adjust as needed)
y_limits <- c(-1, 30)  # Example limits for the y-axis (adjust as needed)


# PC1 Plot
pc1_WL_current_infection <- ggpredict(model_6, terms = c("PC1", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of maximum weight loss") +
    theme_minimal(base_size = 12) +  # Global base font size
    labs(color = "Infection group", fill = "Infection group") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
        legend.position = "none",
        title = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold"),  # X-axis title size
        axis.title.y = element_text(size = 20, face = "bold"),  # Y-axis title size
        axis.text.x = element_text(size = 20),  # X-axis tick size
        axis.text.y = element_text(size = 20),  # Y-axis tick size
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12)
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)

pc1_WL_current_infection

ggsave(paste0(an_fi, "/pc1_WL_current_infection.pdf"), pc1_WL_current_infection, 
       width = 8, height = 6, dpi = 300)  # Use high dpi for publication

# PC2 Plot
pc2_WL_current_infection <- ggpredict(model_6, terms = c("PC2 [-8:5]", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of maximum weight loss") +
    theme_minimal(base_size = 12) +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(color = "Infection group") +
    theme(
        legend.position = "none",
        title = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)

pc2_WL_current_infection

ggsave(paste0(an_fi, "/pc2_WL_current_infection.pdf"), 
       pc2_WL_current_infection, width = 8, height = 6, dpi = 300)  # High dpi





combined_plot <- (
    eimeria_weight /
         biplot /  
        (pc1_WL_current_infection + pc2_WL_current_infection)
) +
    plot_annotation(
        title = "Figure 1",
        theme = theme(
            plot.title = element_text(hjust = 0, size = 16, face = "bold")
        )
    ) +
    plot_layout(guides = "collect") +  # Collect all legends
    plot_annotation(
        tag_levels = "A"  # Automatically add panel labels A, B, C, etc.
    )


ggsave(
    filename = paste0(panels_fi, "/Figure1_combined.pdf"),
    plot = combined_plot,
    width = 12, height = 16, dpi = 300  # Adjusted dimensions
)

 










