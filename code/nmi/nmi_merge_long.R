# --------------------------------------------------
# Script: nmi_merge_long.R
# Purpose: Join cleaned lab (Challenge) and field datasets
# Input: Challenge (lab), field (field)
# Output: merge_prior_imputation.csv
# --------------------------------------------------

# 1. Check overlap between column names -----------------------------

length(intersect(colnames(Challenge), colnames(field)))  # Expected ~37
length(outersect(colnames(Challenge), colnames(field)))  # Expected ~148

# 2. Harmonize factor levels before merge ---------------------------

Challenge$MC.Eimeria <- as.factor(Challenge$MC.Eimeria)

# 3. Merge datasets by shared columns -------------------------------

hm <- full_join(Challenge, field, 
                by = intersect(colnames(field), colnames(Challenge)))

# 4. Clean up after merge -------------------------------------------

# Remove duplicated columns ending in _N (from lab or field redundancy)
hm <- hm %>% dplyr::select(-ends_with("_N"))

# Ensure consistent Mouse_ID format
hm$Mouse_ID <- str_replace(hm$Mouse_ID, "_", "")

# 5. Save merged dataset --------------------------------------------

write.csv(hm, file = paste0(danalysis, "/merge_prior_imputation.csv"), 
          row.names = FALSE)
