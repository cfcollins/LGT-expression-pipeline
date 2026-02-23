### sequence analysis 

## load libraries
library(dplyr)
library(ggplot2)

## load data
pseudogenes <- read.csv("005-degenerating_categorisation_2_truncation/01-data/pseudogene_LGT_vs_native.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession")

# combine with categories
categories <- read.csv("4. How do LGTs compare to Natives and Donors/5. Categories/sum_all_categories.csv") %>% 
  select(LGT.Cluster, LGT_native_accession, Category) %>% 
  filter(Category != "All homologs silenced") %>% 
  mutate(Overall_category = ifelse(Category %in% c("Lower", "Silenced"), "Degenerating", "Retained")) %>% 
  arrange(LGT.Cluster, LGT_native_accession) %>% 
  as_tibble()

  
##### STATS

### LGT vs Native

# make pseudogene a factor
pseudogenes$Pseudogene <- factor(pseudogenes$Pseudogene, levels = c("No", "Yes"))
# run a GLMM
glmm_model_LN2 <- glmer(Pseudogene ~ Type + (1 | LGT.Cluster/Type), data = pseudogenes, family = binomial)
summary(glmm_model_LN2)
sum(resid(glmm_model_LN2, type = "pearson")^2) / 465
plot(resid(glmm_model_LN2))
qqnorm(resid(glmm_model_LN2))
exp(fixef(glmm_model_LN2))

# visualise these results
pseudogenes$predicted_prob <- predict(glmm_model_LN2, type = "response")

# Plot predicted probabilities by Overall_category
ggplot(pseudogenes, aes(x = Type, y = predicted_prob, fill = Type)) +
  geom_boxplot() +
  labs(title = "Predicted Probability of Being a Pseudogene by Overall Category", 
       y = "Predicted Probability", 
       x = "Overall Category") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),     # Axis titles
    axis.text = element_text(size = 12),      # Axis labels
    strip.text = element_text(size = 16),     # Facet labels
    plot.title = element_text(size = 18),     # Plot title
    legend.title = element_text(size = 14),   # Legend title (key)
    legend.text = element_text(size = 12)        # Plot title
  )

ggsave("7. Regulation of LGTs and Natives - Sequence Variation/LGT_vs_native_predicted_probs.png")
### Degenerating vs Retained 
LGTs_categories <- categories %>% 
  rename(Accession = "LGT_native_accession") %>% 
  left_join(pseudogenes, b = c("LGT.Cluster", "Accession")) %>% 
  filter(Type == "LGT") %>% 
  select(LGT.Cluster, Accession, Gene, Overall_category, Pseudogene) 
write.csv(LGTs_categories, "7. Regulation of LGTs and Natives - Sequence Variation/sum_all_categories_pseudogenes.csv", row.names = FALSE)

# make pseudogene a factor
LGTs_categories$Pseudogene <- factor(LGTs_categories$Pseudogene, levels = c("No", "Yes"))

# chi squared test
category_data <- LGTs_categories %>%
  group_by(Overall_category) %>% 
  summarise(Yes_Count = sum(Pseudogene == "Yes", na.rm = TRUE),
            No_Count = sum(Pseudogene == "No", na.rm = TRUE))

chi_squared_result <- chisq.test(category_data[, c("Yes_Count", "No_Count")])
print(chi_squared_result)


# Run a GLMM
glmm_model1 <- glmer(
  Pseudogene ~ Overall_category + (1 | LGT.Cluster),
  data = LGTs_categories,
  family = binomial
)
summary(glmm_model1)

glmm_model2 <- glmer(
  Pseudogene ~ Overall_category + (1 | Accession),
  data = LGTs_categories,
  family = binomial
)
summary(glmm_model2)
table(LGTs_categories$Accession)

anova(glmm_model1, glmm_model2, test = "Chisq")


LGTs_categories$predicted_prob <- predict(glmm_model1, type = "response")

# Plot predicted probabilities by Overall_category
ggplot(LGTs_categories, aes(x = Overall_category, y = predicted_prob, fill = Overall_category)) +
  geom_boxplot() +
  labs(title = "Predicted Probability of Being a Pseudogene by Overall Category", 
       y = "Predicted Probability", 
       x = "Overall Category")


### compared to logCount 
counts <- read.csv("1. Data Preparation/3. Generation of LGT-Native-Donor counts/COUNT_LISTS/LGT_natives_counts.csv") %>% 
  as_tibble()

pseudogenes

combined <- left_join(pseudogenes, counts, by = c("LGT.Cluster", "Accession", "Type", "Gene")) %>% 
  select(LGT.Cluster, Accession, Gene, Type, Pseudogene, logCount)

model <- lmer(
  logCount ~ Pseudogene + (1 | LGT.Cluster/Type),
  data = combined
)
summary(model)

combined %>% 
  ggplot(aes(x = Pseudogene, y = logCount, fill = Pseudogene)) +
  geom_boxplot() +
  facet_wrap(~Type)

lgts <- combined %>% 
  filter(Type == "LGT")
natives <- combined %>% 
  filter(Type == "Native")
lgt_model <- lmer(
  logCount ~ Pseudogene + (1 | LGT.Cluster),
  data = lgts
)
summary(lgt_model)

lgt_model <- lmer(
  logCount ~ Pseudogene + (1 | LGT.Cluster),
  data = natives
)
summary(lgt_model)

lgts %>% 
  ggplot(aes(x = Pseudogene, y = logCount, fill = Pseudogene)) +
  geom_boxplot()

###
pseudogene_counts <- pseudogenes %>% 
group_by(Type, Accession) %>% 
  count(Pseudogene) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
two_colours = c("#8DA0CB", "#E78AC3")
ggplot(pseudogene_counts, aes(x = Type, y = pct * 100, fill = Pseudogene)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%")), 
            position = position_stack(vjust = 0.5),
            size = 5) +
  facet_wrap(~Accession, nrow = 1) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = two_colours) +
  labs(y = "Category Percentage", x = "Accession")
ggsave("figure3b.svg", dpi = 1200, width = 10, height = 5)
