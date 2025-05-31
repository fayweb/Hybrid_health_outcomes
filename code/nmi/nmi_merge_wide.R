# --------------------------------------------------
# Script: nmi_merge_wide.R
# Purpose: Create a unified dataset of experimental mice 
#          (field + selected lab infection groups only)
# Input: gmf (field genes), gml (lab genes), lab, field
# Output: merge_experimental_data.csv
# --------------------------------------------------

# 1. Filter: keep only mice present in gene expression tables -------

field <- field %>% 
    filter(Mouse_ID %in% gmf$Mouse_ID)

lab <- lab %>%
    filter(Mouse_ID %in% gml$Mouse_ID)

# 2. Select lab mice used in experiments -----------------------------

# Primary infection group
death_prim <- lab %>% 
    filter(death == "primary")

# Mice that died during the challenge infection
death_challenge <- lab %>%
    filter(death == "challenge", infection == "challenge")

# Combine selected lab mice
lab <- rbind(death_prim, death_challenge)

# Harmonize Eimeria infection coding
lab$MC.Eimeria <- as.factor(lab$MC.Eimeria)

# Clean workspace
rm(death_prim, death_challenge)

# 3. Merge with cleaned field data -----------------------------------

hm <- rbind(lab, field)

# 4. Save unified dataset for normalization/imputation --------------

write.csv(hm, paste0(danalysis, "/intermediate/merge_experimental_data.csv"),
          row.names = FALSE)

          