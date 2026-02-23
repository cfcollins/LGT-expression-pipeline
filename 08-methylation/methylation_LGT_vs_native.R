#### ANALYSING METHYLATION DATA ####
# LGT vs vertically inherited genes

## Load Libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(svglite)

## load data
lgt_native_meth <- read.csv("007-methylation/01-Data/lgt_native_meth_combined.csv") %>% 
  as_tibble()

## visualise
custom_labels <- c("CDS" = "Coding Sequence", 
                   "full_gene" = "Full Gene")

two_colours = c("#e31f1fff", "#37c8abff")
ggplot(lgt_native_meth, aes(x = Type, y = proportion)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 3, aes(color = Type)) +
  facet_wrap(~Meth_coverage, labeller = labeller(Meth_coverage = custom_labels)) +
  labs(x = "Overall Category", y = "Proportion") +
  scale_colour_manual(values = two_colours) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),     # Axis titles
    axis.text = element_text(size = 14),      # Axis labels
    strip.text = element_text(size = 18),     # Facet labels
    plot.title = element_text(size = 20),     # Plot title
    legend.title = element_text(size = 16),   # Legend title (key)
    legend.text = element_text(size = 14)        # Plot title
  )

ggsave("LGT_vs_Native_meth.svg")

# mixed-effects models

## LGT vs Native
full_gene <- lgt_native_meth %>% 
  filter(Meth_coverage == "full_gene")
cds <- lgt_native_meth %>% 
  filter(Meth_coverage == "CDS")
upstream <- lgt_native_meth %>% 
  filter(Meth_coverage == "1000bp_upstream")
downstream <- lgt_native_meth %>% 
  filter(Meth_coverage == "1000bp_downstream")

# full gene
full_gene_model <- lmer(
  proportion ~ Type + (1 | LGT.Cluster),
  data = full_gene
)
summary(full_gene_model)
anova(full_gene_model)

# cds
cds_model <- lmer(
  proportion ~ Type + (1 | LGT.Cluster/Type),
  data = cds
)
summary(cds_model)
anova(cds_model)

# upstream
upstream_model <- lmer(
  proportion ~ Type + (1 | LGT.Cluster/Type),
  data = upstream
)
summary(upstream_model)
anova(upstream_model)

# downstream
model <- lmer(
  proportion ~ Type + (1 | LGT.Cluster/Type),
  data = downstream
)
summary(model)
anova(model)

## difference
diff <- data %>% 
  group_by(Meth_coverage, LGT.Cluster, Accession, Type) %>% 
  summarise(proportion = mean(proportion)) %>% 
  group_by(Meth_coverage, LGT.Cluster, Accession) %>% 
  summarise(
    difference_in_proportion = proportion[Type == "LGT"] - proportion[Type == "Native"]
  )

full_gene <- diff %>% 
  filter(Meth_coverage == "full_gene")
cds <- diff %>% 
  filter(Meth_coverage == "CDS")
upstream <- diff %>% 
  filter(Meth_coverage == "1000bp_upstream")
downstream <- diff %>% 
  filter(Meth_coverage == "1000bp_downstream")

ggplot(downstream, aes(x = difference_in_proportion)) +
  geom_boxplot()

t.test(downstream$difference_in_proportion, alternative = "greater")

mean(diff$difference_in_proportion)
median(diff$difference_in_proportion)
