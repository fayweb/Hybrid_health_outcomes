#7.2: Aplication of random forest on field samples
# requires hm and random forest model
# Creates Field with updated new variable of predicted weight loss for each mouse
# removed one mouse due to missing genotyping

# filter for the field mice
Field <- hm %>%
    filter(origin == "Field")# %>%
 #   drop_na(HI) # remove any nas in the hybrid index


# select the gene columns
gene <-  Field %>%
    ungroup() %>%
    dplyr::select(c(Mouse_ID, all_of(Genes_v)))

# data frame with only the genes
genes <- gene %>%
    dplyr::select(-Mouse_ID)


# load predicting weight loss model
weight_loss_predict <- readRDS("code/models/WL_predict_gene.rds")

set.seed(540)

#The predict() function in R is used to predict the values 
#based on the input data.
predicted_WL <- predict(weight_loss_predict, genes)

# assign test.data to a new object, so that we can make changes
result_field <- genes

#add the new variable of predictions to the result object
result_field <- cbind(result_field, as.data.frame(predicted_WL))

# add it to the field data 
Field <- cbind(Field,  as.data.frame(predicted_WL))

rm(gene,genes, result_field, weight_loss_predict)

