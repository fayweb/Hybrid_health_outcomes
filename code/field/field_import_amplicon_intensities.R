
# add missing infection intensities
ii <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/CEWE_FECES_infection_intensities")

ii$Mouse_ID <- gsub(pattern = "AA_", replacement = "AA", x = ii$Mouse_ID)


field <- field %>%
    left_join(ii, by = intersect(colnames(field), colnames(ii)))


rm(ii)
