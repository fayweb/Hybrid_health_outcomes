# Requires the file: master_hm.R to run
# File located in: GitHub -> Hybrid_health_outcomes -> code 
# Requires code of master file from line 1 - 163
# Producing a csv file with data on 336 mice
immune_data <- hm %>%
    filter(origin == "Field")

write.csv(x = immune_data, file = "output/tables/field_immune_data.csv", 
          row.names = FALSE)


## these variables are relevant to your analysis
# They have been process in the master script: master_hm.R
#### vectors for selecting genes for analysis
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
               "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
               "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
               "TICAM1", "TNF")