#### categorising putatively stable

## load libraries
library(dplyr)
library(ggplot2)
library(FSA)
library(broom)
library(tidyr)

## load data
put_stab <- read.csv("expanded_deg_or_stab_TOP_GENE_LEAF_SHEATH.csv") %>% 
  as_tibble() %>% 
  filter(deg_or_stab == "Putatively stable") %>% 
  select(-Expanded_Category, -deg_or_stab)

counts <- read.csv("LGT_native_donor_counts_filtered.csv") %>% 
  as_tibble() %>% # process with method
  filter(Sample.organ != "Root") %>% 
  group_by(LGT.Cluster, Type, LGT_native_accession, Accession, Sample.organ, Sample) %>% 
  summarise(logCount = max(logCount)) %>% 
  ungroup() %>% 
  select(-Sample.organ, -Sample, -Accession)

# combine
combined <- left_join(put_stab, counts, by = c("LGT.Cluster", "LGT_native_accession")) 

# compare groups
kruskal_results <- combined %>% 
  group_by(LGT.Cluster, LGT_native_accession) %>% 
  do(tidy(kruskal.test(logCount ~ Type,
                       data = .,)))
p_values <- kruskal_results$p.value
corrected_p_values <- p.adjust(p_values, method = "BH")
kruskal_results$corrected_p_value <- corrected_p_values
kruskal_results <- kruskal_results %>% 
  mutate(Kruskal_wallis_test = ifelse(corrected_p_value < 0.05, "Significant", "Non-significant")) 

## dunn test
# filter to only include the groups where Kruskal-Wallis is significant
sig <- kruskal_results %>% 
  filter(Kruskal_wallis_test == "Significant")
counts_kruskal_sig <- semi_join(combined, sig, by = c("LGT.Cluster", "LGT_native_accession"))

# then perform dunn's test
dunn <- lapply(split(counts_kruskal_sig, paste(counts_kruskal_sig$LGT.Cluster, counts_kruskal_sig$LGT_native_accession)), 
               function(x) dunnTest(data=x, logCount ~ Type, method = 'bh'))
# extract cols and make matrix
p.adj <- as.data.frame(sapply(dunn, function(x) x$res$P.adj))
lengths_of_elements <- lengths(sapply(dunn, function(x) x$res$P.adj))
comp <- as.data.frame(sapply(dunn, function(x) x$res$Comparison))
dunn_results <- cbind(
  pivot_longer(p.adj, cols = everything(), names_to = 'LGT', values_to = 'p.adj'), 
  pivot_longer(comp, cols = everything(), names_to = 'LGT2', values_to = 'comparison')
) %>% 
  select(-LGT2) %>% 
  as_tibble() %>% 
  separate(LGT, into = c("LGT.Cluster", "LGT_native_accession"), sep = " ", remove = TRUE) %>% 
  mutate(Dunn_test = ifelse(p.adj < 0.05, "Significant", "Non-significant")) %>% 
  arrange(LGT.Cluster, LGT_native_accession)

## combine together
kruskal_nonsig <- kruskal_results %>% 
  filter(Kruskal_wallis_test == "Non-significant") %>% 
  select(LGT.Cluster, LGT_native_accession, Kruskal_wallis_test) %>% 
  mutate(Expression = "Expressed")

dunn_results_edit <- dunn_results %>% 
  select(!p.adj) %>% 
  pivot_wider(id_cols = c(LGT.Cluster, LGT_native_accession), 
              names_from = comparison, 
              values_from = Dunn_test) %>% 
  rename(Dunn_donor_LGT = `Donor - LGT`,
         Dunn_donor_native = `Donor - Native`,
         Dunn_LGT_native = `LGT - Native`) %>% 
  mutate(Expression = "Expressed",
         Kruskal_wallis_test = NA)

all <- rbind(kruskal_nonsig, dunn_results_edit) %>% 
  arrange(LGT.Cluster, LGT_native_accession)  

all_categories <- all %>% 
  mutate(Category = case_when(
    Expression == "Silenced" ~ "Silenced",
    Expression == "All homologs silenced" ~ "All homologs silenced",
    Kruskal_wallis_test == "Non-significant" ~ "Same",
    Dunn_donor_native == "Significant" & Dunn_LGT_native == "Significant" & Dunn_donor_LGT == "Non-significant" ~ "Donor",
    Dunn_donor_native == "Significant" & Dunn_LGT_native == "Non-significant" & Dunn_donor_LGT == "Significant" ~ "Native",
    Dunn_donor_native == "Significant" & Dunn_LGT_native == "Non-significant" & Dunn_donor_LGT == "Non-significant" ~ "Same",
    Dunn_donor_native == "Non-significant" & Dunn_LGT_native == "Non-significant" & Dunn_donor_LGT == "Significant" ~ "Same",
    Dunn_donor_native == "Non-significant" & Dunn_LGT_native == "Significant" & Dunn_donor_LGT == "Non-significant" ~ "Same",
    TRUE ~ NA_character_     
  )) %>%
  mutate(Kruskal_wallis_test  = ifelse(is.na(Kruskal_wallis_test ), "Significant", Kruskal_wallis_test))

# figuring out categories
counts <- read.csv("LGT_native_donor_counts_filtered.csv") %>% 
  group_by(LGT.Cluster, Type, LGT_native_accession, Accession, Sample.organ, Sample) %>% 
  summarise(logCount = sum(logCount))
df <- all_categories %>% 
  filter(is.na(Category))
counts %>% 
  semi_join(df, by = c("LGT.Cluster", "LGT_native_accession")) %>% 
  ggplot(aes(x = Type, y = logCount, fill = Accession)) +
  facet_wrap(~LGT.Cluster, scales = "free_y") +
  geom_boxplot()

# edit it based on the plots above
all_categories %>%
  print(n = 51)
all_categories$Category[1] <- "Intermediate"
all_categories$Category[12] <- "Intermediate"
all_categories$Category[16] <- "Higher"
all_categories$Category[21] <- "Intermediate"
all_categories$Category[26] <- "Intermediate"
all_categories$Category[45] <- "Intermediate"


write.csv(all_categories, "putatively_stable_categories_top_gene_leaf_sheath.csv", row.names = FALSE)

# visualise
summary <- all_categories %>% 
  group_by(LGT_native_accession) %>% 
  count(Category) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)

summary$Category  <- factor(summary$Category , levels = c("Higher", "Intermediate",
                                                          "Donor", "Native", "Same",
                                                          "Lower"))
ggplot(summary, aes(x = LGT_native_accession, y = pct * 100, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  scale_fill_manual(values = c(
    "Higher" = "#A8E6A2",  
    "Intermediate" = "#5CCF7F", 
    "Donor" = "#2BAE98",  
    "Native" = "#1F88C0",
    "Same" = "#1563D3",
    "Lower" = "#F46D43"
  )) +
  theme_minimal()

ggsave("putatively_stable_categories_top_gene_leaf_sheath.png")

## stats tests
overall <- summary %>% 
  group_by(Category) %>% 
  summarise(n = sum(n))
overall
# same the highest category
# Define the observed counts and total counts
observed_counts <- c(1,5,15,9,19)
categories <- c("Higher", "Intermediate", "Donor", "Native", "Same")
total_counts <- sum(observed_counts)

# Check the proportion for the "Same" category
same_count <- observed_counts[categories == "Same"]

# Perform a binomial test
binom_test <- binom.test(same_count, total_counts, p = 1/5, alternative = "greater")
binom_test

# donor the highest category (same removed)
# Define the observed counts and total counts
observed_counts <- c(1,5,15,9)
categories <- c("Higher", "Intermediate", "Donor", "Native")
total_counts <- sum(observed_counts)

# Check the proportion for the "Same" category
same_count <- observed_counts[categories == "Donor"]

# Perform a binomial test
binom_test <- binom.test(same_count, total_counts, p = 1/4, alternative = "greater")
binom_test

