#Normalization


# select genes
genes <- hm[, Genes_v]

genes_matrix <- as.matrix(genes)


genes_matrix <- limma::normalizeQuantiles(genes_matrix)


genes <- as.data.frame(genes_matrix)
