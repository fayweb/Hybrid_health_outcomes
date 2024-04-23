
# select genes
genes <- hm[, c(Genes_v, "PPIB", "GAPDH")]

genes_matrix <- as.matrix(genes)


genes_matrix <- limma::normalizeQuantiles(genes_matrix)


genes <- as.data.frame(genes_matrix)
