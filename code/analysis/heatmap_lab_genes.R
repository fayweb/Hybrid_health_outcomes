# Heatmap
### repeating the heatmap on the now imputed data

# turn the data frame into a matrix and transpose it. We want to have each cell 
# type as a row name 
gene <- t(as.matrix(data.frame(mouse_id,genes)))

# turn the first row into column names
gene %>%
  row_to_names(row_number = 1) -> heatmap_data

heatmap_data <- as.data.frame(heatmap_data)

table(rowSums(is.na(heatmap_data)) == nrow(heatmap_data))


# turn the columns to numeric other wise the heatmap function will not work
heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

# remove columns with only NAs 
heatmap_data <- Filter(function(x)!all(is.na(x)), heatmap_data) 

#remove rows with only Nas
heatmap_data <-  heatmap_data[, colSums(is.na(heatmap_data)) != 
                                nrow(heatmap_data)]



#Prepare the annotation data frame
annotation_df <- as_tibble(lab) %>%
  dplyr::select(c("Mouse_ID",  "WL_max", "current_infection", "MC.Eimeria",
                  "delta_ct_cewe_MminusE")) 

annotation_df <- unique(annotation_df) 

annotation_df <- as.data.frame(annotation_df)
annotation_df$MC.Eimeria <- as.factor(annotation_df$MC.Eimeria)

### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$Mouse_ID


# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(heatmap_data)

#remove the unecessary column
annotation_df <- annotation_df %>% dplyr::select(-Mouse_ID, )


# Define colors for each parasite
parasite_colors <- c("E_falciformis" = "coral2",
                     "E_ferrisi" = "chartreuse4",
                     "uninfected" = "cornflowerblue")
# Define your own colors
#my_colors <- colorRampPalette(c("red", "white", "blue"))(100)

# Generate the heat map
heatmap_eim <-
    pheatmap(heatmap_data, annotation_col = annotation_df,# color = my_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_colors = list(current_infection = parasite_colors, 
                                  MC.Eimeria = c("TRUE" = "red",
                                                 "FALSE" = "blue"))) # use annotation_colors



ggsave(filename = paste0(an_fi, "/heatmap_lab_genes.jpeg"), 
       plot = heatmap_eim, 
       width = 8, height = 4, dpi =1000)

rm(annotation_df, gene, genes, heatmap_data, heatmap_eim, mouse_id)
