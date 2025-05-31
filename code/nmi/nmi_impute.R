# --------------------------------------------------
# Script: nmi_impute.R
# Purpose: Impute missing gene expression values (MICE)
# Input: genes (normalized), hm (merged metadata)
# Output: imputed_clean_data.csv, diagnostic plots
# --------------------------------------------------

# ***********************************************************
# Part 1: Assess missing data --------------------------------
# ***********************************************************

# Count missing values per gene
sapply(genes, function(x) sum(is.na(x)))

# Optional: check lab and field subsets
sapply(gl, function(x) sum(is.na(x)))
sapply(gf, function(x) sum(is.na(x)))

# ***********************************************************
# Part 2: Visualize missing data -----------------------------
# ***********************************************************

# Create MICE structure to explore patterns
init <- mice(genes, maxit = 0)
meth <- init$method  # just in case you want to customize methods

# Aggregation plot
jpeg(paste0(fi, "/aggregation_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
aggr(genes, col = c("navyblue", "red"), numbers = TRUE, sortVars = TRUE,
     labels = names(genes), cex.axis = .7, gap = 3,
     ylab = c("Histogram of missing data", "Pattern"))
dev.off()

# Margin plots for selected genes
jpeg(paste0(fi, "/margin_plots.jpeg"), width = 16, height = 12, units = "in", res = 300)
par(mfrow = c(2, 2))
marginplot(genes[, c("IFNy", "IRGM1")])
marginplot(genes[, c("IL.6", "IRGM1")])
marginplot(genes[, c("TICAM1", "IRGM1")])
marginplot(genes[, c("MUC5AC", "IRGM1")])
par(mfrow = c(1, 1))
dev.off()

# Remove gene with too much missing data
genes <- genes[, !(names(genes) %in% "IL.10")]

# ***********************************************************
# Part 3: Imputation using MICE ------------------------------
# ***********************************************************

igf <- mice(genes, m = 5, seed = 500)  # default method unless customized via 'meth'
summary(igf)

# ***********************************************************
# Part 4: Plot imputation diagnostics ------------------------
# ***********************************************************

# Convergence diagnostics
#jpeg(paste0(fi, "/igf_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#plot(igf)
#dev.off()

# XY plots of imputation results
#jpeg(paste0(fi, "/xy_1_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#xyplot(igf, IRGM1 ~ IFNy + CXCR3 + IL.6 + IL.13 + IL1RN + CASP1 + CXCL9 + IDO1 + MPO,
#       pch = 18, cex = 1)
#dev.off()

#jpeg(paste0(fi, "/xy_2_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#xyplot(igf, IRGM1 ~ MUC2 + MUC5AC + MYD88 + NCR1 + PRF1 + RETNLB + SOCS1 + TICAM1 + TNF,
#       pch = 18, cex = 1)
#dev.off()

# Strip plot (distribution of imputations per variable)
#jpeg(paste0(fi, "/stirrplot_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#stripplot(igf, pch = c(20, 21), cex = 1.2)
#dev.off()

# Boxplot of standard deviation across imputations
#jpeg(paste0(fi, "/bw_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#bwplot(igf)
#dev.off()

# Density plot of distributions
#jpeg(paste0(fi, "/density_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
#densityplot(igf)
#dev.off()

# ***********************************************************
# Part 5: Merge imputed data back into metadata --------------
# ***********************************************************

# Extract complete imputed dataset (first of 5)
complete_genes <- complete(igf, 1)

# Multiply by -1 so that higher values = higher expression
complete_genes[] <- -1 * complete_genes

# Attach Mouse_ID and merge back into metadata
result <- data.frame(Mouse_ID = hm$Mouse_ID, complete_genes)

hm_imp <- hm %>%
    dplyr::select(-c(all_of(Genes_v), GAPDH, PPIB)) %>%
    left_join(result, by = "Mouse_ID")

# Optional: check differences
outersect(colnames(hm_imp), colnames(hm))

# Save cleaned, imputed dataset
write.csv(hm_imp, paste0(danal_final, "/imputed_clean_data.csv"), row.names = FALSE)

# Overwrite main object
hm <- hm_imp

# ***********************************************************
# Part 6: Clean up environment -------------------------------
# ***********************************************************

rm(complete_genes, gmf, gf, gl, gml, hm_imp, init, igf, result, 
   genes_matrix, genes, field, lab)
