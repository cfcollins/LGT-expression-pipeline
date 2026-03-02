### dN/dS for truncated genes

## load libraries
library(dplyr)
library(ggplot2)

## load data
dn_ds <- read.csv("009-dn_ds/01-Data/LGT_natives_dn_ds.csv") %>% 
  as_tibble()
truncated <- read.csv("005-degenerating_categorisation_2_truncation/01-data/pseudogene_LGT_vs_native.csv") 

combined <- full_join(dn_ds, truncated, by = c("LGT.Cluster", "Type")) %>% 
  filter(!is.na(Pseudogene),
         !is.na(BM_Significance))

## visualise
summary <- combined %>% 
  group_by(LGT_native_accession, Type, Pseudogene) %>% 
  count(Classification) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
ggplot(summary, aes(x = Pseudogene, y = pct, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  facet_grid(LGT_native_accession ~ Type) +
  theme_minimal()

ggsave("009-dn_ds/04-Degenerating_vs_put_stable/01-sum-all/dn_ds_deg_vs_stab_sum_all.png")

edit <- combined %>% 
  select(LGT.Cluster, Classification, LGT_native_accession, deg_or_stab) %>% 
  group_by(LGT.Cluster) %>% 
  summarise(
    Classification = first(Classification),
    deg_or_stab = case_when(
      all(deg_or_stab == "Degenerating") ~ "Degenerating",
      all(deg_or_stab == "Putatively stable") ~ "Putatively stable",
      TRUE ~ "Mixed"
    )
  )
summary_edit <- edit %>% 
  group_by(deg_or_stab) %>% 
  count(Classification) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
ggplot(summary_edit, aes(x = deg_or_stab, y = pct, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  theme_minimal()

ggsave("009-dn_ds/04-Degenerating_vs_put_stable/01-sum-all/dn_ds_deg_vs_stab_accession_combined_sum_all.png")
