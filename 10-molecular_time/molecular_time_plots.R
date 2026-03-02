library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)

# load data
age_by_cluster <- read.csv("010-molecular_time/ages_by_LGT_cluster.csv") %>% 
  as_tibble() %>% 
  rename(age = "Oldest") %>% 
  filter(!is.na(age))
age_by_accession <- read.csv("010-molecular_time/Copy of ages_by_accession_summary.csv") %>% 
  as_tibble() %>% 
  rename(age = "Median")
low_vs_not <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/06-LGTs_less_than_native_and_donor/01-sum_all/silenced_lower_vs_not_sum_all.csv")
deg_stab <- read.csv("005-degenerating_categorisation_2_truncation/05-combine_with_expression_degenerating/01-sum_all/expanded_deg_or_stab_SUM_ALL.csv") %>% 
  rename(Accession = "LGT_native_accession")
truncated <- read.csv("005-degenerating_categorisation_2_truncation/01-data/pseudogene_LGT_vs_native.csv")
counts <- read.csv("001-Data/002-Counts/LGTs_only_counts.csv")

## lower / silenced vs comparable
low_com <- full_join(age_by_cluster, low_vs_not, by = c("LGT.Cluster")) %>% 
  filter(!is.na(LGT_sig_less_don),
         !is.na(age)) %>% 
  select(LGT.Cluster, age, LGT_native_accession, Category)
ggplot(low_com, aes(x = Category, y = age, fill = Category)) +
  geom_boxplot() +
  facet_wrap(~LGT_native_accession) +
  theme_minimal()

## degenerating vs putatively stable
deg_stab1 <- full_join(age_by_cluster, deg_stab, by = c("LGT.Cluster")) %>% 
  filter(!is.na(deg_or_stab),
         !is.na(age)) %>% 
  select(LGT.Cluster, age, Accession, deg_or_stab)
ggplot(deg_stab1, aes(x = deg_or_stab, y = age, fill = deg_or_stab)) +
  geom_boxplot() +
  facet_wrap(~Accession) +
  theme_minimal()

deg_stab1 <- full_join(age_by_accession, deg_stab, by = c("LGT.Cluster", "Accession")) %>% 
  filter(!is.na(deg_or_stab),
         !is.na(age)) %>% 
  select(LGT.Cluster, age, Accession, deg_or_stab)
ggplot(deg_stab1, aes(x = deg_or_stab, y = age, fill = deg_or_stab)) +
  geom_boxplot() +
  facet_wrap(~Accession) +
  theme_minimal()
ggplot(deg_stab1, aes(x = deg_or_stab, y = age, fill = deg_or_stab)) +
  geom_boxplot() +
  theme_minimal()

## truncation vs non-truncated
trunc <- full_join(age_by_cluster, truncated, by = c("LGT.Cluster")) %>% 
  filter(!is.na(Pseudogene),
         !is.na(age)) %>% 
  select(LGT.Cluster, age, LGT_native_accession, Pseudogene)
ggplot(trunc, aes(x = Pseudogene, y = age, fill = Pseudogene)) +
  geom_boxplot() +
  facet_wrap(~LGT_native_accession) +
  theme_minimal()

## counts vs age
counts_age <- full_join(age_by_cluster, counts, by = "LGT.Cluster") %>% 
  filter(!is.na(age),
         !is.na(logCount))
ggplot(counts_age, aes(x = age, y = Count)) +
  geom_point()
means <- counts_age %>% 
  group_by(LGT.Cluster, age, Accession, Gene) %>% 
  summarise(mean_logCount = mean(logCount),
            se = sd(logCount) / sqrt(n()))
ggplot(means, aes(x = age, y = mean_logCount)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_logCount - se, ymax = mean_logCount + se), width = 0.2) +
  theme_minimal() +
  labs(x = "Age", y = "Mean Log Count") 
model <- lmer(Count ~ age + (1|LGT.Cluster/Gene), data = counts_age)
summary(model)
library(glmmTMB)

model <- glmmTMB(Count ~ age + (1 | LGT.Cluster), 
                 data = counts_age, 
                 family = tweedie(link = "log"))
summary(model)

# Get residuals
residuals_data <- data.frame(
  residuals = residuals(model),
  fitted = fitted(model),  # Predicted values
  age = counts_age$age
)

# Plot residuals vs fitted values
ggplot(residuals_data, aes(x = log(fitted+1), y = log(residuals+1))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Fitted Values", y = "Residuals") +
  ggtitle("Residuals vs Fitted Values")

# Alternatively, you can plot residuals against age:
ggplot(residuals_data, aes(x = age, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Age", y = "Residuals") +
  ggtitle("Residuals vs Age")
# Histogram of residuals
ggplot(residuals_data, aes(x = residuals)) +
  geom_histogram(binwidth = 0.5) +
  theme_minimal() +
  labs(x = "Residuals", y = "Frequency") +
  ggtitle("Histogram of Residuals")

# Q-Q plot
qqnorm(residuals_data$residuals)
qqline(residuals_data$residuals, col = "red")


### ancient rare vs recent ones
age_by_cluster

ggplot(age_by_cluster, aes(x = age)) +
  geom_histogram(bins = 7)

# Test if the median molecular age is significantly lower than expected (e.g., 50)
wilcox.test(age_by_cluster$age, mu = median(age_by_cluster$age), alternative = "less")  

library(moments)
skewness(age_by_cluster$age)

df <- full_join(age_by_cluster, low_vs_not, by = "LGT.Cluster") %>% 
  filter(!is.na(Category),
         !is.na(age))
ggplot(df, aes(x = age, fill = Category)) +
  geom_histogram() +
  facet_wrap(~Category)

mean(age_by_cluster$age)
median(age_by_cluster$age)
summary(age_by_cluster$age)

library(tseries)
jarque.bera.test(age_by_cluster$age)

prop_low <- mean(age_by_cluster$age < median(age_by_cluster$age))
prop_high <- mean(age_by_cluster$age > median(age_by_cluster$age))
prop_low / prop_high  # If >1, more young genes than old

age_by_cluster %>% 
  summarise(a = sum(age <= 1),
            b = sum(age > 1 & age <= 2),
            c = sum(age > 2 & age <= 3),
            d = sum(age > 3 & age <= 4),
            e = sum(age > 4 & age <= 5),
            f = sum(age > 5 & age <= 6),
            g = sum(age > 6))
