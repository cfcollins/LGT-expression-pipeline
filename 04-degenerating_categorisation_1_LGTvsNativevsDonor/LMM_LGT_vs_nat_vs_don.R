### linear models for LGT vs Native vs Donor
### Cat Collins - 02-2025

## libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

## data
data <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/02-making_usable_spreadsheet/LGT_native_donor_counts_edit.csv") %>% 
  as_tibble() %>% 
  filter(Type %in% c("LGT", "Native", "Donor")) %>% 
  group_by(LGT.Cluster, LGT_native_accession, Type, Accession, Sample.organ, Sample) %>% 
  summarise(logCount = sum(logCount))  %>% 
  mutate(LGT_or_homolog = ifelse(Type == "LGT", "LGT", "homologs"))

## visualise
means <- data %>% 
  group_by(LGT.Cluster, LGT_native_accession, Type, LGT_or_homolog, Accession) %>% 
  summarise(mean_logCount = mean(logCount))
ultimate_mean <- data %>% 
  group_by(Type) %>% 
  summarise(mean_logCount = mean(logCount))
ggplot(means, aes(x = Type, y = mean_logCount, fill = LGT_or_homolog)) +
  geom_violin(alpha = 0.5, outlier.shape = NA) +
  geom_point(data = ultimate_mean, aes(x = Type, y = mean_logCount, fill = Type), pch = 4, size = 5) +
  scale_fill_manual(values = c("homologs" = "#005F5F", "LGT" = "#5E2BFF")) +
  theme_minimal() +
  labs(title = "Expression Differences Between LGT and Homologs",
       x = "Gene Type", y = "logCount") +
  theme(legend.position = "none")
ggsave("004-degenerating_categorisation_1_LGTvsNativevsDonor/03-LMM/LGT_native_donor.svg", width = 8, height = 6)

ggplot(means, aes(x = Type, y = mean_logCount, fill = LGT_or_homolog)) +
  geom_violin(
    aes(fill = LGT_or_homolog, fill = after_scale(colorspace::lighten(fill, .5))), 
    size = 1.2, bw = .4
  ) +
  geom_point(data = ultimate_mean, aes(x = Type, y = mean_logCount, fill = Type), 
             pch = 4, size = 5) +
  geom_boxplot(
    fill = "white", size = 1.2, width = .05, outlier.size = 1
  ) +
  scale_fill_manual(values = c("homologs" = "#005F5F", "LGT" = "#5E2BFF")) +
  theme_minimal() +
  labs(title = "Expression Differences Between LGT and Homologs",
       x = "Gene Type", y = "logCount") +
  theme(legend.position = "none")

## linear models 

# homologs (native + donor) vs LGT
model <- lmer(logCount ~ LGT_or_homolog + (1 | LGT.Cluster / LGT_native_accession), data = data)
summary(model)

# LGT vs Native and LGT vs Donor
# first change the ordering so the linear model compares both native and donor
# to LGT
data$Type <- factor(data$Type, levels = c("LGT", "Native", "Donor"))
# linear mixed effects model accounting for clustering by LGT.Cluster and Accession
model <- lmer(logCount ~ Type + (1 | LGT.Cluster / Accession), data = data)
summary(model)
anova(model)

data$Type <- factor(data$Type, levels = c("Native", "LGT", "Donor"))
# linear mixed effects model accounting for clustering by LGT.Cluster and Accession
model <- lmer(logCount ~ Type + (1 | LGT.Cluster / Accession), data = data)
summary(model)
anova(model)

######
data
categories
combined <- full_join(data,categories, by = c("LGT.Cluster", "LGT_native_accession"))
deg <- combined %>% 
  filter(deg_or_stab == "Degenerating")
stab <- combined %>% 
  filter(deg_or_stab == "Putatively stable")

deg$Type <- factor(deg$Type, levels = c("Donor", "LGT", "Native"))
model <- lmer(logCount ~ Type + (1 | LGT.Cluster / Accession), data = deg)
summary(model)
anova(model)

ggplot(deg, aes(x = Type, y = logCount, fill = LGT_or_homolog)) +
  geom_violin(
    aes(fill = LGT_or_homolog, fill = after_scale(colorspace::lighten(fill, .5))), 
    size = 1.2, bw = .4
  ) +
  geom_point(data = ultimate_mean, aes(x = Type, y = mean_logCount, fill = Type), 
             pch = 4, size = 5) +
  geom_boxplot(
    fill = "white", size = 1.2, width = .05, outlier.size = 1
  ) +
  scale_fill_manual(values = c("homologs" = "#005F5F", "LGT" = "#5E2BFF")) +
  theme_minimal() +
  labs(title = "Expression Differences Between LGT and Homologs",
       x = "Gene Type", y = "logCount") +
  theme(legend.position = "none")
stab$Type <- factor(stab$Type, levels = c("Native", "LGT", "Donor"))
model <- lmer(logCount ~ Type + (1 | LGT.Cluster / Accession), data = stab)
summary(model)
anova(model)
