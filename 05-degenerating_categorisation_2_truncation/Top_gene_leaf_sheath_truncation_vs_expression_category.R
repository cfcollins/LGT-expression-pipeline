## truncation analysis - comparing the current degenerating vs stable categories

## load libraries
library(dplyr)
library(ggplot2)

## load data
pseudogenes <- read.csv("v2_pseudogene_LGT_vs_native.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession") %>% 
  filter(Type == "LGT") %>% 
  rename(LGT_native_accession = "Accession")

categories <- read.csv("silenced_lower_vs_not_top_gene_leaf_sheath.csv") %>% 
  as_tibble() %>% 
  select(LGT.Cluster, LGT_native_accession, Category)

data <- full_join(pseudogenes, categories, by = c("LGT.Cluster", "LGT_native_accession")) %>% 
  filter(!is.na(Category),
         !is.na(Pseudogene))

# visualisation
summary <- data %>% 
  group_by(LGT_native_accession , Category) %>% 
  count(Pseudogene) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n) 
ggplot(summary, aes(x = Category, y = pct * 100, fill = Pseudogene)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  facet_grid(~LGT_native_accession, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = .8))

ggsave("truncation_vs_expression_categories_top_gene_leaf_sheath.png")

# make pseudogene a factor
data$Category <- factor(data$Category, levels = c("Degenerating", "Putatively stable"))
# Run a GLMM
glmm_model1 <- glmer(
  Category ~ Pseudogene + (1 | LGT.Cluster),
  data = data,
  family = binomial
)
summary(glmm_model1)

