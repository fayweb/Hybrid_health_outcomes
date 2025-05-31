# Load final cleaned field dataset
field <- read.csv(paste0(dfield_final, "/field_cleaned_data.csv"))

# Gene expression correlation heatmap
gmf <- field[, c("Mouse_ID", Genes_v, "GAPDH")] 
gmf <- gmf %>% select_if(~!all(is.na(.)))
gmf <- gmf[!apply(is.na(gmf[-1]), 1, all), ]
gf <- gmf[, -1]

gene_correlation <- as.matrix(cor(gf, use = "pairwise.complete.obs"))
p.mat <- cor.mtest(gene_correlation)

jpeg(paste0(an_fi, "/cor_genes_field.jpeg"), width = 8, height = 6, units = "in", res = 300)
corrplot(gene_correlation,
         method = "circle",
         tl.col = "black", tl.srt = 45,
         col = brewer.pal(n = 8, name = "RdYlBu"),
         order = "hclust",
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         addCoef.col = "black",
         number.cex = 0.5,
         title = "Field")
dev.off()
rm(gene_correlation, p.mat)

