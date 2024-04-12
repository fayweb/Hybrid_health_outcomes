#select the genes and lab mice
lab <- lab %>%
    dplyr::filter(Position == "mLN") %>%
    group_by(Mouse_ID, infection) %>%
    filter(dpi == max_dpi)


# select the genes, the mouse identifier and the house keeping gene
gml <- lab[, c("Mouse_ID", Genes_v, "PPIB")] 

gml <- unique(gml)

# Remove columns with only NA values
gml <- gml %>% select_if(~!all(is.na(.)))

#remove rows with only nas 
gml <- gml[!apply(is.na(gml[-1]), 1, all), ]

gl <- gml[,-1]

# draw correlation between the genes
gene_correlation <- as.matrix(cor(gl, 
                                  use="pairwise.complete.obs"))

# matrix of the p-value of the correlatio
p.mat <- cor.mtest(gene_correlation)

jpeg(paste0(an_fi, "/cor_genes_lab.jpeg"), 
     width = 8, height = 6, units = "in", res = 300)

    corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust", #hclust reordering
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         addCoef.col = 'black',
         number.cex=0.5, 
         title = "Lab") 
#Add significance level to the correlogram
#remove the values that are insignificant

dev.off()


