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
                        scale = 0.9, alpha = 0.8, point_shape = 21, 
                        point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
                 position = position_nudge(x = 0.2)) +
    theme_minimal() +
    scale_fill_manual(values = color_mapping) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain")  -> eimeria_weight

eimeria_weight

italics_y(eimeria_weight, labels)

ggsave(filename = paste0(d_fi,"/eimeria_weight_combined.jpeg"),
       plot = eimeria_weight, 
       width = 8, height = 6, dpi = 1000)

###############################################      Plot 1B
lab <- hm %>%
    filter(origin == "Lab")

# select the laboratory genes
genes <- lab[ ,colnames(lab) %in% Genes_v]

# PCA
##pca on the complete data set
res.pca <- PCA(genes)

############### biplot
biplot <- fviz_pca_biplot(res.pca, 
                          col.ind = lab$current_infection,
                          pointsize = 2,
                          addEllipses = TRUE, 
                          alpha.ind = 0.9,
                          alpha.var = 0.9, 
                          label = "var",
                          col.var = "black", 
                          repel = TRUE,
                          legend.title = "Infection groups",
                          title = "") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    scale_shape_manual(values = c(15, 16, 17), labels = labels) +
    labs(color = "Infection groups", shape = "Infection groups") +
    theme(legend.text = element_markdown()) 
   # stat_ellipse(geom = "polygon", alpha = 0.6) # Set ellipse alpha


biplot


ggsave(filename = paste0(an_fi, "/biplot.jpeg"), plot = biplot, 
       width = 12, height = 6, dpi = 600)


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

model_3 <- lm(WL_max ~ PC1 + PC2 , data = lab)

summary(model_3)


coef_plot_PC1PC2 <- plot_coefs(model_3)
coef_plot_PC1PC2

coef_plot_PC1PC2 <-
    coef_plot_PC1PC2 +
    theme_minimal(base_size = 14) + # Use a minimal theme with larger font sizes
    theme(
        axis.title = element_text(face = "bold"), # Bold axis titles
        axis.text = element_text(size = 12), # Larger axis text
        legend.text = element_text(size = 12), # Larger legend text
        plot.title = element_text(hjust = 0.5, face = "bold"), # Center the title
       # panel.grid.major = element_line(color = "gray80"), # Subtle grid lines
        panel.grid.minor = element_blank() # Remove minor grid lines
    ) +
    labs(
        title = element_blank(),
        x = "Estimate", # Customize x-axis label
        y = "Principal Components" # Customize y-axis label
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) + # Highlight zero line
    geom_point(size = 3, color = "darkblue") + # Adjust point size and color
    scale_x_continuous(breaks = seq(-2, 1, by = 0.5)) # Customize x-axis ticks 


coef_plot_PC1PC2



#################### Interaction effects model
model_6 <- lm(WL_max ~ PC1 * current_infection + PC2 *current_infection, 
              data = lab)

summary(model_6)


# Now create the scatter plot using this color mapping
ggpredict(model_6, terms = c("PC1", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    labs(color = "Infection group") +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()) -> pc1_WL_current_infection

pc1_WL_current_infection

ggsave(paste0(an_fi, "/pc1_WL_current_infection.jpeg"), pc1_WL_current_infection, 
       width = 8, height = 6, dpi = 1000)


# Now create the scatter plot using this color mapping
# Now create the scatter plot using this color mapping
ggpredict(model_6, terms = c("PC2", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    scale_color_manual(values = color_mapping, labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    labs(color = "Infection group") +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_markdown()) -> pc2_WL_current_infection

pc2_WL_current_infection

ggsave(paste0(tables, "/pc2_WL_current_infection.jpeg"),
       pc2_WL_current_infection, width = 8, height = 6, dpi = 1000)
