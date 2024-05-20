genes <- hm %>% 
    dplyr::select(all_of(Genes_v))


# Example of performing PCA and aligning subspaces
pca_lab <- prcomp(lab_data)
pca_field <- prcomp(field_data)

# Assuming you align the first k principal components
lab_projected <- pca_lab$x[, 1:k]
field_projected <- pca_field$x[, 1:k]

# Use these projected data for further analysis or model training
