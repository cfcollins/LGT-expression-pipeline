## load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## load data - LGTs with present and expressed native/donor homologs
counts <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/04-removing_silent_homologs/LGT_native_donor_counts_filtered.csv") %>% 
  as_tibble()

# filter out silenced genes
expressed <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/05-silenced_or_expressed/silencedvexpressed_LGT.csv") %>%
  as_tibble %>% 
  filter(Expression == "Expressed") %>% 
  rename(LGT_native_accession = "Accession")
  
expressed_counts <- semi_join(counts, expressed, by = c("LGT.Cluster", "LGT_native_accession")) 

# process with method (e.g. sum or top gene)
counts_method <- expressed_counts %>% 
  group_by(LGT.Cluster, Type, LGT_native_accession, Accession, Sample.organ, Sample) %>% 
  summarise(logCount = sum(logCount)) %>% 
  ungroup()

# Perform stats tests per group
# reduce data
df <- counts_method %>% 
  select(LGT.Cluster, LGT_native_accession, Type, logCount)

results <- df %>%
  mutate(Count = as.numeric(logCount)) %>% 
  group_by(LGT.Cluster, LGT_native_accession) %>%
  nest() %>%
  mutate(
    test_LGT_vs_Native = map(data, ~wilcox.test(.x$logCount[.x$Type == "LGT"], 
                                                .x$logCount[.x$Type == "Native"], 
                                                paired = TRUE, 
                                                alternative = "less")$p.value),
    test_LGT_vs_Donor = map(data, ~wilcox.test(.x$logCount[.x$Type == "LGT"], 
                                               .x$logCount[.x$Type == "Donor"], 
                                               paired = FALSE, 
                                               alternative = "less")$p.value)
  ) %>%
  unnest(cols = c(test_LGT_vs_Native, test_LGT_vs_Donor))

# Adjust p-values for multiple testing 
test_LGT_vs_Native  <- results$test_LGT_vs_Native 
p_LGT_vs_Native_adj <- p.adjust(test_LGT_vs_Native, method = "BH")
results$p_LGT_vs_Native_adj <- p_LGT_vs_Native_adj

test_LGT_vs_Donor   <- results$test_LGT_vs_Donor  
p_LGT_vs_Donor_adj <- p.adjust(test_LGT_vs_Donor, method = "BH")
results$p_LGT_vs_Donor_adj <- p_LGT_vs_Donor_adj

# View the results
print(results)

results_edit <- results %>% 
  select(-test_LGT_vs_Native, -test_LGT_vs_Donor, -data) %>% 
  mutate(LGT_sig_less_nat = p_LGT_vs_Native_adj < 0.05,
         LGT_sig_less_don = p_LGT_vs_Donor_adj < 0.05) %>% 
  mutate(LGT_lower_than_homologs = ifelse(
    LGT_sig_less_nat == TRUE & LGT_sig_less_don == TRUE, 
    "Sig_lower", "Not_sig_lower"))

write.csv(results_edit, "004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/LGTs_less_than_native_and_donor_or_not-SUM-ALL.csv", row.names = FALSE)

####
# combine with silenced data
expressed_silenced <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/05-silenced_or_expressed/silencedvexpressed_LGT.csv") %>%
  as_tibble %>% 
  rename(LGT_native_accession = "Accession")

# make it equivalent (e.g. combination of paralogs within accessions)
reduced <- expressed_silenced %>%
  group_by(LGT.Cluster, LGT_native_accession) %>% 
  summarise(Expression = ifelse(any(Expression == "Expressed"), "Expressed", "Silenced"))

# combine with the rest of the data
full <- full_join(reduced, results_edit, by = c("LGT.Cluster", "LGT_native_accession")) %>% 
  mutate(LGT_lower_than_homologs = ifelse(is.na(LGT_lower_than_homologs), "Silenced", LGT_lower_than_homologs)) %>% 
  group_by(LGT_native_accession) %>% 
  mutate(Category = ifelse(LGT_lower_than_homologs %in% c("Silenced", "Sig_lower"), "Degenerating", "Putatively stable"))
write.csv(full, "004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/silenced_lower_vs_not_sum_all.csv", row.names = FALSE)  

# summary stats
summary <- full %>% 
  group_by(LGT_native_accession) %>% 
  summarise(
    n = n(),
    degenerating_n = sum(Category == "Degenerating"),
    degenerating_percen = (degenerating_n / n)*100
  )

write.csv(summary, "004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/summary_stats_LGT_nat_don_SUM_ALL.csv", row.names = FALSE)

# plot
summary_2 <- full %>% 
  group_by(LGT_native_accession) %>% 
  count(Category) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)

summary_2$Category  <- factor(summary_2$Category , levels = c("Putatively stable", "Degenerating"))
summary_2$LGT_native_accession  <- factor(summary_2$LGT_native_accession , levels = c("AUS1", "ZAM", "L04", "KWT", "MRL"))
ggplot(summary_2, aes(x = LGT_native_accession, y = pct * 100, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  theme_minimal()

ggsave("004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/LGT_less_than_homologs_or_not_sum_all.svg")
