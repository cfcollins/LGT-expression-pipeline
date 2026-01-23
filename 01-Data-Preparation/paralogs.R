library(FSA)
library(dplyr)
library(ggplot2)

# accessing count data
counts <- read.csv("01-Data/02-Counts/LGT_natives_counts.csv") %>% 
  as_tibble() %>% 
  filter(Type == "LGT" | Type == "Native") %>% 
  select(!c(Donor, Count, Species, Replicate))

# finding clusters/accessions/types where paralogs exist
all_paralogs <- counts %>% 
  group_by(LGT.Cluster, Accession, Type) %>%
  summarise(UniqueGenes = n_distinct(Gene)) %>% 
  filter(UniqueGenes > 1)

# labelling all these combos on whether paralogs exist or not
summary <- counts %>% 
  group_by(LGT.Cluster, Accession, Type) %>%
  summarise(UniqueGenes = n_distinct(Gene)) %>% 
  mutate(Paralogs = ifelse(UniqueGenes == 1, "No", "Yes"))

# finding percentage of groupings with paralogs for LGT
LGT <- summary %>% 
  filter(Type == "LGT") %>% 
  group_by(Accession) %>%
  summarise(No_paralogs = sum(Paralogs == "No"), 
            Paralogs = sum(Paralogs == "Yes"), 
            percentage = Paralogs / (Paralogs + No_paralogs)*100) %>% 
  mutate(Type = "LGT")
mean(LGT$percentage)
sd(LGT$percentage)
# finding percentage of groupings with paralogs for Native
Native <- summary %>% 
  filter(Type == "Native") %>% 
  group_by(Accession) %>%
  summarise(No_paralogs = sum(Paralogs == "No"), 
            Paralogs = sum(Paralogs == "Yes"), 
            percentage = Paralogs / (Paralogs + No_paralogs)*100) %>% 
  mutate(Type = "Native")
mean(Native$percentage)
sd(Native$percentage)
# adding dataset together
rbind(LGT, Native) %>% 
  select(Accession, percentage, Type) %>% 
  ggplot(aes(x = Type, y = percentage, fill = Accession)) +
  geom_bar(stat = "identity", position = "dodge")
# quick t test to see if there is a difference
wilcox.test(LGT$percentage, Native$percentage)
mean(Native$percentage)

write.csv(summary, "3. Are LGTs Expressed More Than Natives/Dataframes/paralogs.csv", row.names = FALSE)

# two paralogs
paralogs <- counts %>% 
  group_by(LGT.Cluster, Accession, Type) %>%
  summarise(UniqueGenes = n_distinct(Gene)) %>% 
  filter(UniqueGenes == 2) 

paralogs1 <- semi_join(counts, paralogs, by = c("LGT.Cluster", "Accession", "Type"))

paralogs1 <- group_by(paralogs1, LGT.Cluster, Accession, Type)
results <- do(paralogs1, tidy(wilcox.test(.$logCount ~ .$Gene,
                                       alternative = "two.sided")))

p_values <- results$p.value
corrected_p_values <- p.adjust(p_values, method = "hochberg")
results$corrected_p_value <- corrected_p_values
results <- results %>% 
  mutate(Significance = ifelse(is.na(corrected_p_value), "Non-Significant", 
                               ifelse(corrected_p_value < 0.05, "Significant", "Non-Significant"))) %>% 
  print(n=187)
sig <- results %>% 
  filter(Significance == "Significant")
nonsig <- results %>% 
  filter(Significance == "Non-Significant")

ggplot(counts %>% filter(LGT.Cluster == "LGT-001" & Accession == "L04"), aes(x = Type, y = logCount, fill = Gene)) +
  geom_boxplot()

# three or more paralogs 
paralogs_multiple <- counts %>% 
  group_by(LGT.Cluster, Accession, Type) %>%
  summarise(UniqueGenes = n_distinct(Gene)) %>% 
  filter(UniqueGenes > 2)

paralogs_multiple1 <- semi_join(counts, paralogs_multiple, by = c("LGT.Cluster", "Accession", "Type"))

kruskal_results <- paralogs_multiple1 %>% 
  group_by(LGT.Cluster, Accession, Type) %>% 
  do(tidy(kruskal.test(logCount ~ Gene,
                       data = .,)))
p_values <- kruskal_results$p.value
corrected_p_values <- p.adjust(p_values, method = "hochberg")
kruskal_results$corrected_p_value <- corrected_p_values
kruskal_results <- kruskal_results %>% 
  mutate(Kruskal_wallis_test = ifelse(corrected_p_value < 0.05, "Significant", "Non-significant")) 

## dunn test
# filter to only include the groups where Kruskal-Wallis is significant
sig <- kruskal_results %>% 
  filter(Kruskal_wallis_test == "Significant")
counts_kruskal_sig <- semi_join(paralogs_multiple1, sig, by = c("LGT.Cluster", "Accession"))

# then perform dunn's test
dunn <- lapply(
  split(counts_kruskal_sig, paste(counts_kruskal_sig$LGT.Cluster, 
                                               counts_kruskal_sig$Accession, 
                                               counts_kruskal_sig$Type)), 
               function(x) dunnTest(data=x, logCount ~ Gene, method = 'bh'))

lengths_of_elements <- lengths(sapply(dunn, function(x) x$res$P.adj))

#three
three <- dunn[c(1,2,4,6,7,8,12)]
p.adj <- as.data.frame(sapply(three, function(x) x$res$P.adj))
comp <- as.data.frame(sapply(three, function(x) x$res$Comparison))
dunn_results_three <- cbind(
  pivot_longer(p.adj, cols = everything(), names_to = 'LGT', values_to = 'p.adj'), 
  pivot_longer(comp, cols = everything(), names_to = 'LGT2', values_to = 'comparison')) %>% 
  select(-LGT2) %>% 
  as_tibble() %>% 
  separate(LGT, into = c("LGT.Cluster", "Accession", "Type"), sep = " ", remove = TRUE) %>% 
  mutate(Dunn_test = ifelse(p.adj < 0.05, "Significant", "Non-significant")) %>% 
  arrange(LGT.Cluster, Accession)
#six
six <- dunn[c(3,9, 11)]
p.adj <- as.data.frame(sapply(six, function(x) x$res$P.adj))
comp <- as.data.frame(sapply(six, function(x) x$res$Comparison))
dunn_results_six <- cbind(
  pivot_longer(p.adj, cols = everything(), names_to = 'LGT', values_to = 'p.adj'), 
  pivot_longer(comp, cols = everything(), names_to = 'LGT2', values_to = 'comparison')) %>% 
  select(-LGT2) %>% 
  as_tibble() %>% 
  separate(LGT, into = c("LGT.Cluster", "Accession", "Type"), sep = " ", remove = TRUE) %>% 
  mutate(Dunn_test = ifelse(p.adj < 0.05, "Significant", "Non-significant")) %>% 
  arrange(LGT.Cluster, Accession)
#ten
ten <- dunn[c(5,10)]
p.adj <- as.data.frame(sapply(ten, function(x) x$res$P.adj))
comp <- as.data.frame(sapply(ten, function(x) x$res$Comparison))
dunn_results_ten <- cbind(
  pivot_longer(p.adj, cols = everything(), names_to = 'LGT', values_to = 'p.adj'), 
  pivot_longer(comp, cols = everything(), names_to = 'LGT2', values_to = 'comparison')) %>% 
  select(-LGT2) %>% 
  as_tibble() %>% 
  separate(LGT, into = c("LGT.Cluster", "Accession", "Type"), sep = " ", remove = TRUE) %>% 
  mutate(Dunn_test = ifelse(p.adj < 0.05, "Significant", "Non-significant")) %>% 
  arrange(LGT.Cluster, Accession) 

all_multiples <- rbind(dunn_results_three, dunn_results_six, dunn_results_ten)


# combine
results 
kruskal_results %>% 
  filter(Kruskal_wallis_test == "Non-significant")
all_multiples %>% 
  print(n = 59)


## checking the effect of paralogs on method
# checking to see how paralogs affect this
paralogs <- read.csv("3. Are LGTs Expressed More Than Natives/Dataframes/paralogs.csv") %>% 
  as_tibble()
yes <- paralogs %>% 
  filter(Paralogs == "Yes")
no <- paralogs %>% 
  filter(Paralogs == "No")

counts <- read.csv("1. Data Preparation/3. Generation of LGT-Native-Donor counts/COUNT_LISTS/LGT_natives_counts.csv") %>% 
  filter(Type == "LGT" | Type == "Native") %>% 
  as_tibble()

sum_all <- counts %>% 
  group_by(LGT.Cluster, Type, Sample, Accession) %>% 
  summarise(logCount = sum(logCount)) %>% 
  mutate(Method = "Sum All")
top_gene_all <- counts %>%  
  group_by(LGT.Cluster, Type, Sample, Accession) %>% 
  summarise(logCount = max(logCount)) %>% 
  mutate(Method = "Top Gene All")
sum_leaf_sheath <- counts %>%  
  filter(Sample.organ != "Root") %>% 
  group_by(LGT.Cluster, Type, Sample, Accession) %>% 
  summarise(logCount = sum(logCount)) %>% 
  mutate(Method = "Sum Leaf Sheath")
top_gene_leaf_sheath <- counts %>%  
  filter(Sample.organ != "Root") %>% 
  group_by(LGT.Cluster, Type, Sample, Accession) %>% 
  summarise(logCount = max(logCount)) %>% 
  mutate(Method = "Top Gene Leaf Sheath")

#sum_all
yes_counts <- semi_join(sum_all %>% filter(Type == "LGT"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(sum_all %>% filter(Type == "LGT"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#top_gene_all
yes_counts <- semi_join(top_gene_all %>% filter(Type == "LGT"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(top_gene_all %>% filter(Type == "LGT"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#sum_leaf_sheath
yes_counts <- semi_join(sum_leaf_sheath %>% filter(Type == "LGT"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(sum_leaf_sheath %>% filter(Type == "LGT"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#top_gene_leaf_sheath
yes_counts <- semi_join(top_gene_leaf_sheath %>% filter(Type == "LGT"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(top_gene_leaf_sheath %>% filter(Type == "LGT"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)

#sum_all
yes_counts <- semi_join(sum_all %>% filter(Type == "Native"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(sum_all %>% filter(Type == "Native"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#top_gene_all
yes_counts <- semi_join(top_gene_all %>% filter(Type == "Native"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(top_gene_all %>% filter(Type == "Native"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#sum_leaf_sheath
yes_counts <- semi_join(sum_leaf_sheath %>% filter(Type == "Native"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(sum_leaf_sheath %>% filter(Type == "Native"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)
#top_gene_leaf_sheath
yes_counts <- semi_join(top_gene_leaf_sheath %>% filter(Type == "Native"), yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(top_gene_leaf_sheath %>% filter(Type == "Native"), no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$logCount, no_counts$logCount)

# averaging???
average <- rbind(sum_all, top_gene_all, sum_leaf_sheath, top_gene_leaf_sheath) %>% 
  group_by(LGT.Cluster, Accession, Type, Sample) %>% 
  summarise(mean_logCount = mean(logCount))
yes_counts <- semi_join(average, yes, by = c("LGT.Cluster", "Accession", "Type"))
no_counts <- semi_join(average, no, by = c("LGT.Cluster", "Accession", "Type"))
wilcox.test(yes_counts$mean_logCount, no_counts$mean_logCount)

all <- rbind(sum_all, top_gene_all, sum_leaf_sheath, top_gene_leaf_sheath)
combined <- left_join(all, paralogs, by = c("LGT.Cluster", "Accession", "Type"))

ggplot(combined, aes(x = Method, y = logCount, fill = Paralogs)) +
  geom_boxplot() +
  facet_wrap(~ Type)
