
# select the genes, the mouse identifier and the house keeping gene
gmf <- field[, c("Mouse_ID", Genes_v, "GAPDH")] 

# Remove columns with only NA values
gmf <- gmf %>% select_if(~!all(is.na(.)))

#remove rows with only nas 
gmf <- gmf[!apply(is.na(gmf[-1]), 1, all), ]

gf <- gmf[,-1]


# draw correlation between the genes
gene_correlation <- as.matrix(cor(gf, 
                                  use="pairwise.complete.obs"))

# matrix of the p-value of the correlatio
p.mat <- cor.mtest(gene_correlation)

jpeg(paste0(an_fi, "/cor_genes_field.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)

corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust", #hclust reordering
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         addCoef.col = 'black',
         number.cex=0.5, 
         title = "Field") 
#Add significance level to the correlogram
#remove the values that are insignificant

dev.off()

rm(gene_correlation, p.mat)


