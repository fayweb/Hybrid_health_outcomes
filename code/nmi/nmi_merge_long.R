
# Join wild and lab data 
length(intersect(colnames(Challenge), colnames(field)))

#37 intersecting columns


length(outersect(colnames(Challenge), colnames(field)))

# 148 dissimilar columns

#expected columns:
37 + 153 #190

Challenge$MC.Eimeria <- as.factor(Challenge$MC.Eimeria)

# now join the two data sets
hm <- full_join(Challenge, field, 
                  by = intersect(colnames(field), colnames(Challenge)))

hm <- hm %>%
    dplyr::select(-ends_with("_N"))


hm$Mouse_ID <- str_replace(hm$Mouse_ID, "_", "")


write.csv(x = hm, file = paste0(danalysis, "/merge_prior_imputation.csv"), 
          row.names = FALSE)

