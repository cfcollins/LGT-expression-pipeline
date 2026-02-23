# truncation - truncation effect on logCount

## load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(svglite)
library(glmmTMB)

## load data
pseudogenes <- read.csv("005-degenerating_categorisation_2_truncation/01-data/v2_pseudogene_LGT_vs_native.csv") %>% 
  as_tibble() %>% 
  rename(Accession = "LGT_native_accession") %>% 
  select(LGT.Cluster, Accession, Type, Gene, Pseudogene)

counts <- read.csv("001-Data/002-Counts/LGT_natives_counts.csv") %>% 
  as_tibble() %>% 
  select(LGT.Cluster, Accession, Type, Gene, logCount)

data <- full_join(pseudogenes, counts, by = c("LGT.Cluster", "Accession", "Type", "Gene")) %>% 
  filter(!is.na(logCount),
         !is.na(Pseudogene))
library(colorspace)

# visualise
ggplot(data, aes(x = Pseudogene, y = logCount, color = Pseudogene, fill = Pseudogene)) +
  geom_violin(
    aes(fill = Pseudogene, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, bw = .2
  ) + 
  geom_boxplot(
    fill = "white", size = 1.2, width = .05, outlier.size = 1
  ) +
  facet_wrap(~Type) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "Yes" = "#8B5FBF",  # Light purple for LGT pseudogene
    "No" = "#5E2BFF",   # Dark purple for LGT non-pseudogene
    "Native_Yes" = "#4A8572", # Light teal for Native pseudogene
    "Native_No" = "#005F5F"   # Dark teal for Native non-pseudogene
  )) +
  scale_color_manual(values = c(
    "Yes" = "#8B5FBF",  
    "No" = "#5E2BFF",   
    "Native_Yes" = "#4A8572", 
    "Native_No" = "#005F5F"   
  ))

ggsave("005-degenerating_categorisation_2_truncation/03-truncation_effect_on_logCount/violin_plot.svg", width = 10, height = 5)


means <- data %>% 
  group_by(Type, Pseudogene) %>% 
  summarise(logCount = mean(logCount))

ggplot(data, aes(x = Pseudogene, y = logCount, fill = Pseudogene)) +
  geom_violin() +
  facet_wrap(~Type) +
  theme_minimal()

ggplot(data %>% filter(Type == "LGT"), aes(x = Pseudogene, y = logCount, fill = Pseudogene)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~Accession) +
  theme_minimal()

## stats

model1 <- lmer(logCount ~ Pseudogene + Accession + (1|LGT.Cluster),
              data = data)
summary(model1)
model2 <- lmer(logCount ~ Pseudogene + Accession + (1|Accession),
              data = data)
summary(model2)
model3 <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
               data = data)
summary(model3)
model4 <- lmer(logCount ~ Pseudogene + Accession + (1|LGT.Cluster/Accession/Type),
               data = data)
summary(model4)

anova(model1, model2, test = "Chisq")
anova(model1, model3, test = "Chisq")
anova(model3, model4, test = "Chisq")

# by type
# LGT
LGT <- data %>% 
  filter(Type == "LGT")
LGT_model1 <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster),
                  data = LGT)
summary(LGT_model1)

# native
native <- data %>% 
  filter(Type == "Native")
native_model1 <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster),
                   data = native)
summary(native_model1)


# by accession
# AUS1
AUS1 <- LGT %>% 
  filter(Accession == "AUS1")
AUS1_model <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
                   data = AUS1)
summary(AUS1_model)
# KWT
KWT <- LGT %>% 
  filter(Accession == "KWT")
KWT_model <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
                   data = KWT)
summary(KWT_model)
# L04
L04 <- LGT %>% 
  filter(Accession == "L04")
L04_model <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
                   data = L04)
summary(L04_model)
# MRL
MRL <- LGT %>% 
  filter(Accession == "MRL")
MRL_model <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
                   data = MRL)
summary(MRL_model)
# ZAM
ZAM <- LGT %>% 
  filter(Accession == "ZAM")
ZAM_model <- lmer(logCount ~ Pseudogene + (1|LGT.Cluster/Type),
                   data = ZAM)
summary(ZAM_model)


#####

library(ggplot2)
library(lme4)

# Model residuals
resid_plot <- ggplot(data.frame(resid = residuals(model, type = "pearson")), aes(x = resid)) +
  geom_histogram(bins = 30, fill = "gray", color = "black") +
  theme_minimal() +
  labs(title = "Residual Distribution", x = "Residuals", y = "Frequency")

# QQ plot
qq_plot <- ggplot(data.frame(resid = residuals(model)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "QQ Plot of Residuals")

print(resid_plot)
print(qq_plot)

residual_variance <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
print(residual_variance)
VarCorr(model)



model_nb <- glmmTMB(logCount ~ Pseudogene + (1 | LGT.Cluster), 
                    data = data, family = nbinom2)  # nbinom1 (mean~variance), nbinom2 (variance~mean^2)
summary(model_nb)
anova(model1, model_nb, test = "Chisq")
