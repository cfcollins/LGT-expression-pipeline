## creating Figure 3 

## load libraries
library(ggplot2)
library(dplyr)

## LGT vs Native
lgt_native_meth <- read.csv("007-methylation/01-Data/lgt_native_meth_combined.csv") %>% 
  as_tibble()
full_gene <- lgt_native_meth %>% 
  filter(Meth_coverage == "full_gene")

means <- full_gene %>% 
  group_by(Type) %>% 
  summarise(mean_proportion = mean(proportion))
two_colours = c("#5e2bffff", "#005f5fff")
ggplot(full_gene, aes(x = Type, y = proportion)) +
  geom_boxplot(data = means, aes(x = Type, y = mean_proportion)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 3, aes(color = Type)) +
  scale_colour_manual(values = two_colours) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),     
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 18),     
    plot.title = element_text(size = 20),    
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 14)       
  )
ggsave("007-methylation/LGT_vs_native_full_gene.svg")

## silenced/lower vs everything else
categories <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/silenced_lower_vs_not_sum_all.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession") %>% 
  select(LGT.Cluster, Accession, Category)

full <- full_join(full_gene, categories, by = c("LGT.Cluster", "Accession")) %>% 
  filter(!is.na(Category),
         !is.na(Meth_coverage),
         Type == "LGT")

means <- full %>% 
  group_by(Category) %>% 
  summarise(mean_proportion = mean(proportion))
two_colours = c("#f46d43ff", "#1a9850ff")
ggplot(full, aes(x = Category, y = proportion)) +
  geom_boxplot(data = means, aes(x = Category, y = mean_proportion)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 3, aes(color = Category)) +
  scale_colour_manual(values = two_colours) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),     
    axis.text = element_text(size = 14),      
    strip.text = element_text(size = 18),     
    plot.title = element_text(size = 20),    
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 14)       
  )
ggsave("007-methylation/lower_silenced_vs_not_full_gene.svg")
