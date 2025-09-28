library(ggplot2)
library(ggforce)
library(patchwork)
library(dplyr)

# Define consistent color scheme
color_mapping <- c("Uninfected" = "#4ACAFF",
                   "E. ferrisi" = "#7A0092",
                   "E. falciformis" = "#FF0000")

# Panel A: Your actual PCA data (unchanged)
create_actual_pca <- function(data) {
    data <- data %>%
        mutate(infection_group = case_when(
            grepl("uninfected", current_infection, ignore.case = TRUE) ~ "Uninfected",
            current_infection == "E. ferrisi" ~ "E. ferrisi",
            current_infection == "E. falciformis" ~ "E. falciformis",
            TRUE ~ "Other"
        )) %>%
        filter(infection_group != "Other")
    
    p1 <- ggplot(data, aes(x = PC1, y = PC2, color = infection_group)) +
        geom_point(alpha = 0.6, size = 2) +
        stat_ellipse(level = 0.95, size = 1) +
        scale_color_manual(values = color_mapping) +
        
        # Gene vector arrows
        geom_segment(aes(x = 0, y = 0, xend = -1, yend = 3),
                     arrow = arrow(length = unit(0.25, "cm")),
                     color = "darkgreen", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
        geom_segment(aes(x = 0, y = 0, xend = 3.5, yend = 0.5),
                     arrow = arrow(length = unit(0.25, "cm")),
                     color = "darkred", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
        geom_segment(aes(x = 0, y = 0, xend = 2, yend = -2),
                     arrow = arrow(length = unit(0.25, "cm")),
                     color = "darkorange", size = 0.8, alpha = 0.6, inherit.aes = FALSE) +
        
        # Gene labels
        annotate("text", x = -1.2, y = 3.3, label = "SOCS1\nIRGM1\nMUC2",
                 size = 2.3, color = "darkgreen", fontface = "italic") +
        annotate("text", x = 3.8, y = 0.8, label = "TNF\nIDO1\nCXCL9",
                 size = 2.3, color = "darkred", fontface = "italic") +
        annotate("text", x = 2.2, y = -2.3, label = "MPO\nIL1RN",
                 size = 2.3, color = "darkorange", fontface = "italic") +
        
        geom_point(aes(x = 0, y = 0), inherit.aes = FALSE,
                   size = 2, shape = 21, fill = "white", color = "black") +
        
        labs(x = "PC1 (34.2% variance)",
             y = "PC2 (15.8% variance)",
             title = "A",
             color = "Infection") +
        theme_minimal() +
        theme(panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(face = "bold", size = 12),
              legend.position = "bottom",
              panel.grid.minor = element_blank()) +
        coord_cartesian(xlim = c(-6, 5), ylim = c(-6, 5))
    
    return(p1)
}

# Panel B: Functional zones without parasite labels
# Panel B: Fixed with neutral colors
create_disease_space_zones <- function(data) {
    p2 <- ggplot() +
        xlim(-6, 5) + ylim(-6, 5) +
        
        # Functional zones with neutral colors
        annotate("rect", xmin = -3, xmax = 0.5, ymin = -0.5, ymax = 1.5,
                 fill = "gray80", alpha = 0.3) +  # Homeostatic
        
        annotate("rect", xmin = 0, xmax = 2.5, ymin = 0, ymax = 2,
                 fill = "lightgreen", alpha = 0.3) +  # Resolution
        
        annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -2, ymax = 0,
                 fill = "orange", alpha = 0.2) +  # Early inflammatory
        
        annotate("rect", xmin = 2, xmax = 4, ymin = -2, ymax = 1,
                 fill = "salmon", alpha = 0.3) +  # Peak inflammatory
        
        # Zone labels
        annotate("text", x = -1, y = 0.8, label = "Homeostatic\n(baseline immunity)",
                 size = 3, fontface = "italic", color = "gray50") +
        annotate("text", x = 1.2, y = 1.2, label = "Resolution\n(regulatory response)",
                 size = 3, fontface = "italic", color = "darkgreen") +
        annotate("text", x = 0, y = -1.2, label = "Early inflammation\n(neutrophil-driven)",
                 size = 3, fontface = "italic", color = "darkorange") +
        annotate("text", x = 3, y = -0.5, label = "Peak inflammation\n(Th1-driven)",
                 size = 3, fontface = "italic", color = "darkred") +
        
        labs(x = "PC1 (Inflammatory axis)",
             y = "PC2 (Regulatory axis)",
             title = "B.") +
        theme_minimal() +
        theme(panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(face = "bold", size = 12),
              panel.grid.minor = element_blank())
    
    return(p2)
}

# Panel C: Biologically sensible trajectories
# Panel C: Biologically realistic immune trajectories
# Panel C: Simple single immune trajectory loop
create_infection_cycle <- function() {
    # Create a single smooth loop trajectory
    t <- seq(0, 2*pi, length.out = 100)
    
    # Parametric equations for a tilted elliptical loop
    loop <- data.frame(
        x = 2.5 * cos(t) + 1.2 * cos(t) * sin(t),
        y = 1.8 * sin(t) - 0.5
    )
    
    p3 <- ggplot() +
        xlim(-6, 5) + ylim(-6, 5) +
        
        # The main trajectory loop
        geom_path(data = loop,
                  aes(x = x, y = y),
                  color = "gray40", size = 2, alpha = 0.7) +
        
        # Add directional arrow
        geom_segment(aes(x = 2, y = -2, xend = 2.8, yend = -1.5),
                     arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                     color = "gray40", size = 1.5, alpha = 0.7) +
        
        # Origin point (baseline)
        geom_point(aes(x = 0, y = 0),
                   size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
        
        # Sample points representing day 8 sampling at different phases
        geom_point(data = data.frame(x = c(0.5, 2.5, 3, 1, -1),
                                     y = c(-1, -2, 0, 1.2, 0.5)),
                   aes(x = x, y = y),
                   size = 4, color = "red", alpha = 0.8) +
        
        # Simple labels
        annotate("text", x = 0, y = -0.5, label = "Baseline", size = 3, fontface = "bold") +
        annotate("text", x = 3, y = -2.5, label = "Peak inflammation", size = 3, color = "gray40") +
        annotate("text", x = 2, y = 2, label = "Resolution", size = 3, color = "gray40") +
        
        # Bottom annotation
        annotate("text", x = 0, y = -5,
                 label = "Infections follow a trajectory through immune space\nDay 8 sampling (red) captures different phases",
                 size = 3, fontface = "italic", hjust = 0.5) +
        
        labs(x = "PC1 (Inflammatory axis)",
             y = "PC2 (Regulatory axis)",
             title = "C.") +
        theme_minimal() +
        theme(panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(face = "bold", size = 12),
              panel.grid.minor = element_blank())
    
    return(p3)
}

# Panel D: Clean weight loss landscape
# Panel D: Using actual model predictions
create_weight_loss_landscape <- function(data) {
    data_clean <- data %>%
        filter(!is.na(PC1) & !is.na(PC2) & !is.na(WL_max))
    
    # Create prediction grid
    x <- seq(-6, 5, length.out = 50)
    y <- seq(-6, 5, length.out = 50)
    grid <- expand.grid(PC1 = x, PC2 = y)
    
    # Use Model 1 coefficients (PC only model)
    grid$WeightLoss <- 9.66 +        # intercept
        0.51 * grid$PC1 -   # PC1 effect
        1.27 * grid$PC2     # PC2 effect
    
    # Bound predictions to reasonable range
    grid$WeightLoss <- pmax(0, pmin(25, grid$WeightLoss))
    
    p4 <- ggplot(grid, aes(x = PC1, y = PC2)) +
        geom_contour_filled(aes(z = WeightLoss), bins = 10) +
        geom_contour(aes(z = WeightLoss), color = "gray50", alpha = 0.5, size = 0.3) +
        scale_fill_manual(values = rev(RColorBrewer::brewer.pal(10, "RdYlGn")),
                          name = "Predicted\nWeight Loss (%)") +
        
        # Add text annotations in corners
       # annotate("text", x = 4, y = -5, 
        #         label = "High inflammation\nLow regulation\n= Maximum damage",
           #      size = 2.5, color = "darkred", fontface = "italic") +
        #annotate("text", x = -5, y = 4.5,
         #        label = "Low inflammation\nHigh regulation\n= Minimal damage",
          #       size = 2.5, color = "darkblue", fontface = "italic") +
        
        labs(x = "PC1 (Inflammatory axis)",
             y = "PC2 (Regulatory axis)",
             title = "D.") +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              panel.border = element_rect(fill = NA, color = "black"),
              plot.title = element_text(face = "bold", size = 12))
    
    return(p4)
}

# Create all panels
p1 <- create_actual_pca(lab)
p2 <- create_disease_space_zones(lab)
p3 <- create_infection_cycle()
p4 <- create_weight_loss_landscape(lab)



# Combine into final figure
final_figure <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
       # title = "Immune Trajectories and Disease Outcomes",
        #subtitle = "Cross-sectional sampling captures different phases of dynamic immune responses",
        theme = theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 11, face = "italic")
        )
    )

final_figure

# Save outputs
ggsave("immune_trajectory_figure.png", final_figure, width = 12, height = 10, dpi = 300)
ggsave("immune_trajectory_figure.pdf", final_figure, width = 12, height = 10)

