# --------------------------------------------------
# Script: nmi_normalize.R
# Purpose: Apply quantile normalization to gene expression data
# Input: hm (merged lab + field dataset)
# Output: genes (normalized expression data)
# --------------------------------------------------

# 1. Subset gene expression columns (plus housekeeping genes) ----
genes <- hm[, c(Genes_v, "PPIB", "GAPDH")]

# 2. Convert to matrix for normalization --------------------------
genes_matrix <- as.matrix(genes)

# 3. Apply quantile normalization (limma package) -----------------
genes_matrix <- limma::normalizeQuantiles(genes_matrix)

# 4. Convert back to dataframe for downstream use -----------------
genes <- as.data.frame(genes_matrix)

# Note: Save output here if needed (e.g. as intermediate CSV)
# write.csv(genes, paste0(danalysis, "/intermediate/genes_normalized.csv"), row.names = FALSE)
