## molecular time by accession

# libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyr)
library(stringr)

# data
age_by_cluster <- read.csv("010-molecular_time/ages_by_LGT_cluster.csv") %>% 
  as_tibble() %>% 
  rename(age = "Oldest") %>% 
  filter(!is.na(age))
deg_stab <- read.csv("005-degenerating_categorisation_2_truncation/05-combine_with_expression_degenerating/01-sum_all/expanded_deg_or_stab_SUM_ALL.csv") %>% 
  rename(Accession = "LGT_native_accession")
counts <- read.csv("001-Data/002-Counts/LGTs_only_counts.csv")
counts_age <- full_join(age_by_cluster, counts, by = "LGT.Cluster") %>% 
  filter(!is.na(age),
         !is.na(logCount))
age_plus_cat <- full_join(age_by_cluster, deg_stab, by = "LGT.Cluster") %>% 
  filter(!is.na(deg_or_stab),
         !is.na(age)) %>% 
  select(LGT.Cluster, age, Accession, deg_or_stab)

# split up accessions
AUS1_counts_age <- counts_age %>% 
  filter(Accession == "AUS1")
KWT_counts_age <- counts_age %>% 
  filter(Accession == "KWT")
L04_counts_age <- counts_age %>% 
  filter(Accession == "L04")
ZAM_counts_age <- counts_age %>% 
  filter(Accession == "ZAM")
MRL_counts_age <- counts_age %>% 
  filter(Accession == "MRL")

AUS1_cat_age <- age_plus_cat %>% 
  filter(Accession == "AUS1")
KWT_cat_age <- age_plus_cat %>% 
  filter(Accession == "KWT")
L04_cat_age <- age_plus_cat %>% 
  filter(Accession == "L04")
ZAM_cat_age <- age_plus_cat %>% 
  filter(Accession == "ZAM")
MRL_cat_age <- age_plus_cat %>% 
  filter(Accession == "MRL")
# plot time vs logCount (altogether)
means <- counts_age %>% 
  group_by(LGT.Cluster, age, Accession, Gene) %>% 
  summarise(mean_logCount = mean(logCount),
            se = sd(logCount) / sqrt(n()))
ggplot(means, aes(x = age, y = mean_logCount)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_logCount - se, ymax = mean_logCount + se), width = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal() +
  labs(x = "Age", y = "Mean Log Count") 
# plot time vs logCount (by accession)
ggplot(means, aes(x = age, y = mean_logCount)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_logCount - se, ymax = mean_logCount + se), width = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal() +
  labs(x = "Age", y = "Mean Log Count")+
  facet_wrap(~Accession, nrow = 2) +
  theme(
    text = element_text(size = 20),            # increase all text size
    strip.text = element_text(size = 20),      # facet label size
    axis.title = element_text(size = 20),      # axis title size
    axis.text = element_text(size = 20)        # axis tick label size
  )

ggsave("time_vs_logCount_by_accession.png", width = 17, height = 10)

# plot time vs category (altogether)
ggplot(age_plus_cat, aes(x = deg_or_stab, y = age, fill = deg_or_stab)) +
  geom_boxplot() 

# plot time vs category (by accession)
ggplot(age_plus_cat, aes(x = deg_or_stab, y = age, fill = deg_or_stab)) +
  geom_boxplot() +
  geom_jitter(width = .2) +
  facet_wrap(~Accession, nrow =1) +
  scale_fill_manual(values = c(
    "Degenerating" = "#f46d43ff",
    "Putatively stable" = "#1a9850ff"
  )) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),            # increase all text size
    strip.text = element_text(size = 20),      # facet label size
    axis.title = element_text(size = 20),      # axis title size
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),       # axis tick label size
    axis.text.x = element_blank()  # rotate x-axis labels
  )

ggsave("time_vs_category_by_accession.svg", width = 25, height = 4)
# linear model time vs logCount (altogether)
mod1 <- lm(logCount ~ age, data = counts_age)
summary(mod1)

# linear model time vs logCount (by accession)
AUS1_mod1 <- lm(logCount ~ age, data = AUS1_counts_age)
summary(AUS1_mod1)
KWT_mod1 <- lm(logCount ~ age, data = KWT_counts_age)
summary(KWT_mod1)
L04_mod1 <- lm(logCount ~ age, data = L04_counts_age)
summary(L04_mod1)
ZAM_mod1 <- lm(logCount ~ age, data = ZAM_counts_age)
summary(ZAM_mod1)
MRL_mod1 <- lm(logCount ~ age, data = MRL_counts_age)
summary(MRL_mod1)

# mixed effects model time vs logCount (altogether)
mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = counts_age)
summary(mod2)

# mixed effects model time vs logCount (by accession)
AUS1_mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = AUS1_counts_age)
summary(AUS1_mod2)
KWT_mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = KWT_counts_age)
summary(KWT_mod2)
L04_mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = L04_counts_age)
summary(L04_mod2)
ZAM_mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = ZAM_counts_age)
summary(ZAM_mod2)
MRL_mod2 <- lmer(logCount ~ age + (1|LGT.Cluster/Gene), data = MRL_counts_age)
summary(MRL_mod2)

# glm model time vs category (altogether)
age_plus_cat$deg_or_stab <- factor(age_plus_cat$deg_or_stab)
mod3 <- glm(deg_or_stab ~ age, data = age_plus_cat, family = binomial)
summary(mod3)

# linear model time vs category (by accession)
AUS1_mod3 <- glm(deg_or_stab ~ age, data = AUS1_cat_age, family = binomial)
summary(AUS1_mod3)
KWT_mod3 <- glm(deg_or_stab ~ age, data = KWT_cat_age, family = binomial)
summary(KWT_mod3)
L04_mod3 <- glm(deg_or_stab ~ age, data = L04_cat_age, family = binomial)
summary(L04_mod3)
ZAM_mod3 <- glm(deg_or_stab ~ age, data = ZAM_cat_age, family = binomial)
summary(ZAM_mod3)
MRL_mod3 <- glm(deg_or_stab ~ age, data = MRL_cat_age, family = binomial)
summary(MRL_mod3)

wilcox.test(age ~ deg_or_stab, data = age_plus_cat)
wilcox.test(age ~ deg_or_stab, data = AUS1_cat_age)
wilcox.test(age ~ deg_or_stab, data = KWT_cat_age)
wilcox.test(age ~ deg_or_stab, data = L04_cat_age)
wilcox.test(age ~ deg_or_stab, data = ZAM_cat_age)
wilcox.test(age ~ deg_or_stab, data = MRL_cat_age)


# mixed effects model time vs category (altogether)
mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = age_plus_cat, family = binomial)
summary(mod4)

# mixed effects model time vs category (by accession)
AUS1_mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = AUS1_cat_age, family = binomial)
summary(AUS1_mod4)
KWT_mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = KWT_cat_age, family = binomial)
summary(KWT_mod4)
L04_mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = L04_cat_age, family = binomial)
summary(L04_mod4)
ZAM_mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = ZAM_cat_age, family = binomial)
summary(ZAM_mod4)
MRL_mod4 <- glmer(deg_or_stab ~ age + (1|LGT.Cluster), data = MRL_cat_age, family = binomial)
summary(MRL_mod4)



### age plus logCount plus cat

counts <- counts %>%
  mutate(LGT.Cluster = str_replace(LGT.Cluster, "^LGT-042[a]$", "LGT-042"))

combined <- age_by_cluster %>% 
  full_join(counts, by = "LGT.Cluster") %>% 
  select(LGT.Cluster, Accession, Gene, Sample, logCount, age) %>% 
  full_join(deg_stab, by = c("LGT.Cluster", "Accession")) %>% 
  mutate(deg_or_stab = replace_na(deg_or_stab, "No donor")) %>% 
  select(-Expanded_Category) 
  
means <- combined %>% 
  group_by(LGT.Cluster, age, Accession, deg_or_stab, Gene) %>% 
  summarise(mean_logCount = mean(logCount),
            se = sd(logCount) / sqrt(n()))
means %>% 
  print(n = 200)
means <- means %>%
  mutate(deg_or_stab = factor(deg_or_stab,
                              levels = c("Degenerating", "Putatively stable", "No donor")))
ggplot(means, aes(x = age, y = mean_logCount, colour = deg_or_stab)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_logCount - se, ymax = mean_logCount + se), width = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = c("Degenerating" = "#f46d43ff",
                                 "Putatively stable" = "#1a9850ff",
                                 "No donor" = "darkgrey")) +
  theme_minimal() +
  labs(x = "Age", y = "Mean Log Count") +
  facet_wrap(~Accession, nrow = 1)
ggsave("accessions_by_age_logCount_cat.svg")
