# ***********************************************************
# Part 6: Analysis                           ----
# ***********************************************************
#----------------------------------------------------------*
# 6.1: PCA 
# performing a pca analysis on the fieldoratory immune gene data
# Requires: hm
# Creates: field
# Plots: biplot, pca_variables,
#  contr_PC1, contr_PC2
#----------------------------------------------------------*

field <- hm %>%
    filter(origin == "Field") %>%
    drop_na(species_Eimeria)

# select the fieldoratory genes
genes <- field[ ,colnames(field) %in% Genes_v]

# PCA
##pca on the complete data set
res.pca <- PCA(genes)
biplot <- fviz_pca_biplot(res.pca, 
                          col.ind = field$species_Eimeria,
                          pointsize = 2,
                          addEllipses = TRUE, 
                          fieldel = "var",
                          col.var = "black", 
                          repel = TRUE,
                          legend.title = "Parasite",
                          title = "") #+
   # scale_color_manual(values = color_mapping, fieldels = fieldels) +
    #scale_fill_manual(values = color_mapping, fieldels = fieldels) +
    #scale_shape_manual(values = c(15, 16, 17), fieldels = fieldels) +
    #fields(color = "Infection groups", shape = "Infection groups") +
    #theme(legend.text = element_markdown())


biplot



## How much does each dimension contribute to variance?
jpeg(paste0(an_fi, "/scree_plot_pca_field.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)

fviz_eig(res.pca, addfieldels = TRUE, ylim = c(0, 70), barfill = "seagreen2")

dev.off()


jpeg(paste0(an_fi, "/contributions_dimensions_field.jpeg"), 
     width = 10, height = 5, units = "in", res = 300)

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"),
             repel = TRUE, title = "") 
dev.off()

#fviz_pca_ind(res.pca, col.ind = "cos2", 
#            gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"), 
#           repel = TRUE, title = "")


## Description of the dimensions
## We get a correlation between each variable and the first dimension
dimdesc(res.pca)


# Convert mouse_id to a data frame
mouse <- data.frame(Mouse_ID = field[,1])

# Add the new column pc1 to the mouse_id data frame
mouse$PC1 <- res.pca$ind$coord[, 1]
mouse$PC2 <- res.pca$ind$coord[, 2]  # indexing the second column
mouse$PC3 <-  res.pca$ind$coord[, 3]
mouse$PC4 <-  res.pca$ind$coord[, 4]
mouse$PC5 <-  res.pca$ind$coord[, 5]

field <- field %>% 
    left_join(mouse, by = "Mouse_ID")


## We also need to extract the data for the variable contributions to each of 
## the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source(paste0(c, "/functions.R"))

circ <- circleFun(c(0,0),2,npoints = 500)


#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- as.data.frame(res.pca$var$contrib)
var.contrib.matrix <- data.matrix(var.contrib)
corrplot(var.contrib.matrix, is.corr=FALSE) 

pca_var <- as.data.frame(pca.vars)


### Contributions to the first dimension

# Contributions of variables to PC1
jpeg(paste0(an_fi, "/contributions_pc1_field.jpeg"), 
     width = 6, height = 4, units = "in", res = 300)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "Contribution of immune genes to the first dimension of the PCA", 
             fill =  "seagreen2") 

dev.off()


fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "Contribution of immune genes to the first dimension of the PCA", 
             fill =  "seagreen2") -> contributions_pc1 
# res.pca$var$contrib


### Contributions to the second dimension

## Contributions of variables to PC2
jpeg(paste0(an_fi, "/contributions_pc2_field.jpeg"), 
     width = 6, height = 4, units = "in", res = 300)

fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "Contribution of immune genes to the second dimension of the PCA",
             fill =  "seagreen2") 

dev.off()

fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "Contribution of immune genes to the second dimension of the PCA",
             fill =  "seagreen2") -> contributions_pc2 

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)

# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

#select same rows in the first table
field <- field[row.names(genes), ]

##############################
#########
###########################
vpg <- pca.vars

# Change the first column of the variance contribution of variables to the gene
# names
vpg <- vpg %>%
    dplyr::rename(Variable = vars, PC1 = Dim.1, PC2 = Dim.2)

field$cos2 <- field$PC1^2 + field$PC2^2

# PCA graph of individuals
pca_individuals <-
    ggplot(field, aes(x = PC1, y = PC2, color = current_infection)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 5, alpha = 0.5, color = "black", 
                   shape = 21, aes(fill = current_infection))   +
    scale_fill_manual(values = color_mapping, fieldels = fieldels) +
        fields(x = "PC1 (34.37%)", y = "PC2 (16.03%)",# title = "PCA graph of individuals",
         fill = "Treatment group") +
    theme_minimal() +
    theme(#plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_markdown()) 
#guides(color = guide_legend(override.aes = list(size = 4))) 

pca_individuals

ggsave(filename = paste0(an_fi, "/pca_individuals.jpeg"),
       plot = pca_individuals, 
       width = 6, height = 4, dpi = 300)

###########################
####### PCA graph of variables
# Add cos2 variable to the dataframe
vpg$cos2 <- with(vpg, PC1^2 + PC2^2)


# Define custom gradient colors
gradient_colors <- c("#B00B69", "#A55EA7", "#1D1CC9")

# Define the breaks and fieldels for the color legend
breaks <- c(0, 50, 100, 150)
fieldels_2 <- c("0", "50", "100", "150")

# Plotting the factor map 
pca_variables <-
    ggplot(vpg, aes(x = PC1, y = PC2, color = cos2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3) +
    geom_fieldel_repel(aes(fieldel = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    coord_equal() +
    xfield("PC1 (34.37%)") +
    yfield("PC2 (16.03%") +
    #ggtitle("PCA Plot of Variables") +
    theme_minimal()  +
    scale_color_gradient(low = "blue", high = "orange")+
    fields(color = "Squared distance to origin")

pca_variables

ggsave(filename = paste0(an_fi, "/pca_variables.jpeg"), 
       plot = pca_variables, 
       width = 5, height = 6, dpi = 300)

###############################################
##################
############### biplot
biplot <- fviz_pca_biplot(res.pca, 
                col.ind = field$current_infection,
                pointsize = 2,
                addEllipses = TRUE, 
                fieldel = "var",
                col.var = "black", 
                repel = TRUE,
                legend.title = "Infection groups",
                title = "") +
    scale_color_manual(values = color_mapping, fieldels = fieldels) +
    scale_fill_manual(values = color_mapping, fieldels = fieldels) +
    scale_shape_manual(values = c(15, 16, 17), fieldels = fieldels) +
    fields(color = "Infection groups", shape = "Infection groups") +
    theme(legend.text = element_markdown())


biplot


ggsave(filename = paste0(an_fi, "/biplot.jpeg"), plot = biplot, 
       width = 12, height = 6, dpi = 600)

### Create the panel figure
PCA_panel <-
    (pca_variables | biplot ) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add fieldels (A, B, C, etc.)

# Add a figure title
PCA_panel <- PCA_panel + 
    plot_annotation(title = 'Fig. 3', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))



# Save the panel figure
ggsave(paste0(panels_fi, "/panel_pca.jpeg"), 
       PCA_panel, width = 16, height = 8, dpi = 300)


