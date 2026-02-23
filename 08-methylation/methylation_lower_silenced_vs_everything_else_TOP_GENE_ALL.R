## analysing methylation - lower/silenced vs everything else

## load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)

## load data
meth <- read.csv("lgt_native_meth_combined.csv") %>% 
  as_tibble() %>% 
  select(Meth_coverage, LGT.Cluster, Accession, Type, proportion)
categories <- read.csv("silenced_lower_vs_not_top_gene_all.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession") %>% 
  select(LGT.Cluster, Accession, Category)

full <- full_join(meth, categories, by = c("LGT.Cluster", "Accession")) %>% 
  filter(!is.na(Category),
         !is.na(Meth_coverage),
         Type == "LGT")

## visualise 

# Calculate the mean for each category and meth_coverage
mean_data <- full %>%
  group_by(Meth_coverage, Category) %>%
  summarize(mean_proportion = mean(proportion, na.rm = TRUE), .groups = "drop") %>% 
  filter(Meth_coverage %in% c("CDS", "full_gene"))

# Plot with mean line
ggplot(full %>% filter(Meth_coverage %in% c("CDS", "full_gene")), aes(x = Category, y = proportion)) +
  geom_violin() +
  geom_jitter(aes(), size = 1, width = .1) +
  facet_wrap(~Meth_coverage) +
  geom_boxplot(data = mean_data, aes(x = Category, y = mean_proportion), 
             colour = "red", size = .5, width = .3) +  # Red points for means
  theme_minimal()

ggsave("meth_lower_silenced_everything_else_TOP_GENE_ALL.png")

## stats tests

# overall 
full_gene <- full %>% 
  filter(Meth_coverage == "full_gene")
cds <- full %>% 
  filter(Meth_coverage == "CDS")
upstream <- full %>% 
  filter(Meth_coverage == "1000bp_upstream")
downstream <- full %>% 
  filter(Meth_coverage == "1000bp_downstream")
# full_gene
full_gene_model <- lmer(
  proportion ~ Category + (1 | LGT.Cluster),
  data = full_gene
)
summary(full_gene_model)
# CDS
cds_model <- lmer(
  proportion ~ Category + (1 | LGT.Cluster),
  data = cds
)
summary(cds_model)
# upstream
upstream_model <- lmer(
  proportion ~ Category + (1 | LGT.Cluster),
  data = upstream
)
summary(upstream_model)
# downstream
downstream_model <- lmer(
  proportion ~ Category + (1 | LGT.Cluster),
  data = downstream
)
summary(downstream_model)

