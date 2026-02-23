## truncation analysis - comparing the current degenerating vs stable categories

## load libraries
library(dplyr)
library(ggplot2)

## load data
pseudogenes <- read.csv("005-degenerating_categorisation_2_truncation/01-data/v2_pseudogene_LGT_vs_native.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession") %>% 
  filter(Type == "LGT") %>% 
  rename(LGT_native_accession = "Accession")

categories <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/silenced_lower_vs_not_sum_all.csv") %>% 
  as_tibble() %>% 
  select(LGT.Cluster, LGT_native_accession, Category)

data <- full_join(pseudogenes, categories, by = c("LGT.Cluster", "LGT_native_accession")) %>% 
  filter(!is.na(Category),
         !is.na(Pseudogene))

# combine with current degenerating
combined <- data %>% 
  group_by(LGT.Cluster, LGT_native_accession, Category) %>% 
  summarise(Truncated = ifelse(any(Pseudogene == "No"), "No", "Yes")) %>%  
  mutate(Expanded_Category = case_when(
    Truncated == "Yes" & Category == "Degenerating" ~ "Degenerating - truncated and lower expression",
    Truncated == "Yes" & Category == "Putatively stable" ~ "Degenerating - truncated only",
    Truncated == "No" & Category == "Degenerating" ~ "Degenerating - lower expression only",
    TRUE ~ "Putatively stable"
  )) %>% 
  mutate(deg_or_stab = ifelse(Expanded_Category == "Putatively stable", "Putatively stable", "Degenerating")) %>% 
  select(-Truncated, -Category)

write.csv(combined, "expanded_deg_or_stab_SUM_ALL.csv", row.names = FALSE)

# summary stats
summary <- combined %>% 
  group_by(LGT_native_accession) %>% 
  count(Expanded_Category) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n) %>% 
  ungroup()
write.csv(summary, "SUMMARY_expanded_deg_or_stab_SUM_ALL.csv")

# plot
summary$Expanded_Category  <- factor(summary$Expanded_Category , levels = c("Putatively stable", 
                                                                            "Degenerating - truncated only",
                                                                            "Degenerating - lower expression only", 
                                                                            "Degenerating - truncated and lower expression"))
summary$LGT_native_accession <- factor(summary$LGT_native_accession, levels = c("AUS1", "ZAM", "L04", "KWT", "MRL"))
ggplot(summary, aes(x = LGT_native_accession, y = pct * 100, fill = Expanded_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%")), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  scale_fill_manual(values = c(
    "Degenerating - truncated and lower expression" = "#D73027",  # Dark Red
    "Degenerating - truncated only" = "#FDAE61",  # Medium Red-Orange
    "Degenerating - lower expression only" = "#F46D43",  # Light Red-Orange
    "Putatively stable" = "#1A9850"  # Green
  )) +
  theme_minimal()

ggsave("expanded_deg_or_stab_SUM_ALL.svg")
