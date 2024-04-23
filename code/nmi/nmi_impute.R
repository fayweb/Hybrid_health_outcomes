# ***********************************************************
# Title: Predicting the health outcomes of infections in hybrid mice

# Purpose: Impute missing gene expression data
#
# Authors: Fay Webster
# ID variables:
# ***********************************************************
# Part 1: Analyze nas ----
# ***********************************************************
# how many nas does the data set have?
# select genes
sapply(genes, function(x) sum(is.na(x)))

# in the laboratory data set?
sapply(gl, function(x) sum(is.na(x)))

# in the field data set?
sapply(gf, function(x) sum(is.na(x)))

# ***********************************************************
# Part 2: Analysis of missing data structure. ----
# ***********************************************************
# The frequency distribution !f the missing cases per variable can be obtained 
# as:
init <- mice(genes, maxit = 0)

#we want to impute only the specific variables
meth <- init$method

######################### Quality control plots 
#jpeg(paste0(fi, "/aggregation_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)

# Aggregation plot 
#aggr(genes, col=c('navyblue','red'), 
 #                 numbers=TRUE, sortVars=TRUE, 
  #                labels=names(genes), cex.axis=.7, 
   #               gap=3, ylab=c("Histogram of missing data","Pattern"))

#dev.off()

# Creating the margin plots
jpeg(paste0(fi, "/margin_plots.jpeg"), 
     width = 16, height = 12, units = "in", res = 300)
par(mfrow = c(2,2))
plot1 <- marginplot(genes[,c("IFNy", "IRGM1")])
plot2 <- marginplot(genes[,c("IL.6", "IRGM1")])
plot3 <- marginplot(genes[,c("TICAM1", "IRGM1")])
plot4 <- marginplot(genes[,c("MUC5AC", "IRGM1")])

par(mfrow = c(1,1))
dev.off()


# removing il 10 because of abundance of missing data 
genes <- genes[, !(names(genes) %in% "IL.10")]


# ***********************************************************
# Part 3: Multiple imputation with MICE. ----
# ***********************************************************
igf <- mice(genes, m = 5, seed = 500) # method = meth,

summary(igf)
# ***********************************************************

# Part 4: Imputation plots ----
# ***********************************************************
# plot the iterations
jpeg(paste0(fi, "/igf_plot.jpeg"), width = 8, height = 6, units = "in", res = 300)
plot(igf)
dev.off()

jpeg(paste0(fi, "/xy_1_plot.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)
xyplot(igf, IRGM1 ~ IFNy + CXCR3 + IL.6 + IL.13 + IL1RN + CASP1 + CXCL9,
       IDO1 + MPO, pch = 18, cex = 1)
dev.off()


jpeg(paste0(fi, "/xy_2_plot.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)
xyplot(igf, IRGM1 ~  MUC2 + MUC5AC + MYD88 + NCR1 + PRF1 + 
           RETNLB + SOCS1 + TICAM1 + TNF, pch = 18, cex = 1)
dev.off()


#xyplot(igf, IRGM1 ~ IL.13 + PRF1 + CASP1, pch=18, cex=1)

# stirplot
jpeg(paste0(fi, "/stirrplot_plot.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)
stripplot(igf, pch = c(20,21), cex = 1.2)
dev.off()

# plot of sd after each iteration
jpeg(paste0(fi, "/bw_plot.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)
bwplot(igf)
dev.off()

# density plot after each iteration
jpeg(paste0(fi, "/density_plot.jpeg"), width = 8, height = 6, units = "in", 
     res = 300)
densityplot(igf)
dev.off()


#xyplot(igf, IRGM1 ~ IL.13 + MYD88 + IFNy, pch=18, cex=1)
#xyplot(igf, IRGM1 ~ IL.13 + PRF1 + CASP1, pch=18, cex=1)

# ***********************************************************
# Part 5: Merging imputated data set ----
# ***********************************************************
## igf$imp$IFNy
#Now we can get back the completed dataset using the complete()
complete_genes <- complete(igf, 1)

# Transform all columns by multiplying by -1 (higher expression - more positive)
complete_genes[] <- -1 * complete_genes

# join the imputed genes
Mouse_ID <- hm$Mouse_ID
result <- data.frame(Mouse_ID, complete_genes)

hm_imp <- hm %>%
    dplyr::select(-c(all_of(Genes_v), GAPDH, PPIB)) %>%
    left_join(result, by = "Mouse_ID")

outersect(colnames(hm_imp), colnames(hm))

write.csv(hm_imp, paste0(danal_final, "/imputed_clean_data.csv"), row.names = FALSE)

hm <- hm_imp

# remove unecessary data objects
rm(complete_genes, gmf, gf, gl, gml, hm_imp, init, igf, result, 
   genes_matrix, genes, field, lab)

