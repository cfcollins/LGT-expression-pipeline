# truncation - LGT vs Native

## load libraries
library(dplyr)
library(ggplot2)

## load data
pseudogenes <- read.csv("v2_pseudogene_LGT_vs_native.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession")

## visualise 
summary <- pseudogenes %>% 
  group_by(Accession, Type) %>% 
  count(Pseudogene) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
summary %>%
  filter(Pseudogene == "Yes",
         Accession != "MRL") %>% 
  group_by(Type) %>% 
  summarise(
    mean_pct = mean(pct*100),
    SD_pct = sd(pct*100)
  )
pseudogenes %>% 
  group_by(Type) %>% 
  count(Pseudogene) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)

write.csv(summary, "truncation_LGT_vs_native_by_accession.csv", row.names = FALSE)
summary$Accession  <- factor(summary$Accession , levels = c("AUS1", "ZAM", "L04", "KWT", "MRL"))

# by accession
ggplot(summary, aes(x = Type, y = pct * 100, fill = Pseudogene)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%", " ", n)), 
            position = position_stack(vjust = 0.5),
            size = 2.5) +
  facet_grid(~Accession, scales = "free") +
  theme_minimal()
ggsave("truncation_LGT_vs_native_by_accession.svg")

overall_summary <- pseudogenes %>% 
  group_by(Type) %>% 
  count(Pseudogene) %>% 
  mutate(pct = n / sum(n),               
         ypos = cumsum(n) - 0.5 * n)
# overall
ggplot(overall_summary, aes(x = Type, y = pct * 100, fill = Pseudogene)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(sprintf("%1.1f", pct * 100), "%")), 
            position = position_stack(vjust = 0.5),
            size = 2.5) 
  theme_minimal()
  ggsave("truncation_LGT_vs_native_by_overall.png")
  
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

ggsave("LGT_vs_native_predicted_probs.png")
