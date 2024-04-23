##select same rows in the first table
field <- field %>% 
    filter(Mouse_ID %in% gmf$Mouse_ID)


##select same rows in the first table
lab <- lab %>%
    filter(Mouse_ID %in% gml$Mouse_ID)

### select only the mice with gene expreession data

# mice that diedd in the primary infections

death_prim <- lab %>% 
    filter(death == "primary")

# mice that died in the callenge infections
death_challenge <- lab %>%
    filter(death == "challenge", infection == "challenge")

lab <- rbind(death_prim, death_challenge)
lab$MC.Eimeria <- as.factor(lab$MC.Eimeria)

rm(death_prim, death_challenge)

## merge experimental mice
hm <- rbind(lab, field)


write.csv(hm, paste0(danalysis, "/intermediate/merge_experimental_data.csv"),
                     row.names = FALSE)
          
          