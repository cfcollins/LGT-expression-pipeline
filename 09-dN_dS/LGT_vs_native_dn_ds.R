### dN/dS
### comparing branch and branch site models for LGT and vertically inhertited genes

## load libraries
library(dplyr)
library(ggplot2)

## load data
LGT_dnds <- read.csv("009-dn_ds/01-Data/dN_dS_comparisons.csv") %>% 
  as_tibble() %>% 
  mutate(edit = ifelse(Classification %in% c("NS strong purifying selection", "NS weak purifying selection"), "NS purifying selection", Classification)) %>% 
  select(-Classification) %>% 
  rename(Classification = "edit",
         BM_Significance = "Significance") %>% 
  mutate(Type = "LGT") %>% 
  select(LGT.Cluster, BM_Significance, Classification, Type)

LGT_dnds %>% 
  distinct(Classification)

native_BM <- read.csv("009-dn_ds/01-Data/branch_model_native_stats.csv") %>% 
  as_tibble() %>% 
  select(LGT.Cluster, Significance, M0_omega, BM_M0_omega, BM_native_omega) %>% 
  rename(BM_Significance = "Significance")
native_BSM <- read.csv("009-dn_ds/01-Data/branch_site_model_native_stats.csv") %>% 
  as_tibble() %>% 
  select(LGT.Cluster, Significance) %>% 
  rename(BSM_Significance = "Significance")

native_combined <- full_join(native_BM, native_BSM, by = "LGT.Cluster") %>% 
  mutate(Classification = case_when(
    BM_Significance == "Non-significant" ~ "NS purifying selection",
    BM_Significance == "Significant" & M0_omega > BM_native_omega ~ "Stronger purifying selection",
    BM_Significance == "Significant" & M0_omega < BM_native_omega & BSM_Significance == "Non-significant" ~ "Relaxed constraint - no positive selection sites",
    BM_Significance == "Significant" & M0_omega < BM_native_omega & BSM_Significance == "Significant" ~ "Relaxed constraint - positive selection sites"
  )) %>% 
  mutate(Type = "Native") %>% 
  select(LGT.Cluster, BM_Significance, Classification, Type)

combined <- rbind(LGT_dnds, native_combined)

write.csv(combined, "009-dn_ds/01-Data/LGT_natives_dn_ds.csv", row.names = FALSE)

all_summary <- combined %>% 
  group_by(Type) %>% 
  count(Classification) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
## visualise

# pie chart of branch model significance 
data <- combined %>% 
  group_by(Type) %>% 
  count(BM_Significance)

ggplot(data, aes(x = "", y = n, fill = BM_Significance)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  facet_wrap(~Type) +
  theme_void()

ggsave("009-dn_ds/02-LGT_vs_Native/pie_charts_BM_sig.svg")
# bar chart of those that are BM significant 
summary <- combined %>% 
  filter(BM_Significance == "Significant") %>% 
  group_by(Type) %>% 
  count(Classification) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
summary$Classification  <- factor(summary$Classification, levels = c("Relaxed constraint - positive selection sites", 
                                                                "Relaxed constraint - no positive selection sites",
                                                                "Stronger purifying selection"))
ggplot(summary, aes(x = Classification, y = n, fill = Classification)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap( ~Type) +
  theme_minimal()


ggsave("009-dn_ds/02-LGT_vs_Native/bar_charts_BSM_classification.svg", width = 9, height = 5)


LGT_dnds %>% 
  filter(Classification == "Relaxed constraint - positive selection sites")




####
# figuring out categories
categories <- read.csv("005-degenerating_categorisation_2_truncation/05-combine_with_expression_degenerating/01-sum_all/expanded_deg_or_stab_SUM_ALL.csv")
combined_cat <- full_join(combined, categories) %>% 
  filter(BM_Significance == "Significant",
         Type == "LGT")
combined %>% 
  filter(BM_Significance == "Significant")

combined_cat %>% 
  select(LGT.Cluster, Classification, deg_or_stab) %>% 
  print(n = 49)
summary 
combined_cat %>% 
  group_by(LGT.Cluster) %>% 
  summarise(
    Classification = first(Classification),
    cat_edit = case_when(
      all(deg_or_stab == "Degenerating") ~ "All degenerating",
      all(deg_or_stab == "Putatively stable") ~ "All putatively stable",
      all(is.na(deg_or_stab)) ~ "No donor",
      TRUE ~ "Mixed"
    )
  ) %>% 
  group_by(cat_edit) %>% 
  count(Classification) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
summary$Classification  <- factor(summary$Classification, levels = c("Relaxed constraint - positive selection sites", 
                                                                     "Relaxed constraint - no positive selection sites",
                                                                     "Stronger purifying selection"))
