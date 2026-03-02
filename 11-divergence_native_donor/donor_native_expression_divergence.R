######## comparing expression divergence between native and donor genes
##################

## load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
#####################

## load data
count_data <- read.csv("004-degenerating_categorisation_1_LGTvsNativevsDonor/02-making_usable_spreadsheet/LGT_native_donor_counts_edit.csv") %>% 
  as_tibble() %>% 
  filter(Type %in% c("Native", "Donor", "LGT"))

# process with method (e.g. sum or top gene)
count_data_edit <- count_data %>% 
  group_by(LGT.Cluster, Type, LGT_native_accession, Accession, Sample.organ, Sample) %>% 
  summarise(logCount = sum(logCount)) %>% 
  ungroup()

# combine with categories
categories <- read.csv("005-degenerating_categorisation_2_truncation/05-combine_with_expression_degenerating/01-sum_all/expanded_deg_or_stab_SUM_ALL.csv") 

data <- full_join(count_data_edit, categories, by = c("LGT.Cluster", "LGT_native_accession")) 

########################

## create Table S5
differences <- data %>% 
  group_by(LGT.Cluster, LGT_native_accession, deg_or_stab, Type) %>% 
  summarise(mean_logCount = mean(logCount)) %>% 
  filter(Type != "LGT") %>% 
  pivot_wider(id_cols = c(LGT.Cluster, LGT_native_accession, deg_or_stab),
              names_from = Type,
              values_from = mean_logCount) %>% 
  mutate(mean_logCount_diff = Donor - Native,
         Direction = ifelse(mean_logCount_diff > 0, "Donor proxy gene greater expression", "Vertically inherited gene greater expression"),
         Difference = abs(mean_logCount_diff)) %>% 
  rename(Accession = "LGT_native_accession",
         mean_donor = "Donor",
         mean_native = "Native") %>% 
  select(LGT.Cluster, Accession, deg_or_stab, mean_donor, mean_native, Difference, Direction) %>% 
  filter(!is.na(Direction),
         !is.na(deg_or_stab))
differences 

ggplot(differences, aes(x = deg_or_stab, y = log(Difference+1), fill = deg_or_stab)) +
  geom_boxplot() +
  geom_jitter(width = .1, alpha = .5) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Putatively stable" = "#1a9850ff",
      "Degenerating"  = "#f46d43ff"
    )
  )

ggsave("expression_divergence.svg", width = 8, height = 7)

write.csv(differences, "expression_divergence_differences.csv", row.names = FALSE)

model <- lm(log(Difference+1) ~ deg_or_stab, data = differences)
summary(model)
#################

## figure out mean logCount across put and deg, plus SDs
differences %>% 
  group_by(deg_or_stab) %>% 
  summarise(mean_difference = mean(Difference),
            n = n(),
            SD = sd(Difference))

## wilcoxon test for mean differences between put and deg
put <- differences %>% 
  filter(deg_or_stab == "Putatively stable")
deg <- differences %>% 
  filter(deg_or_stab == "Degenerating")
wilcox.test(put$Difference, deg$Difference)

## LMM to test interaction between gene status and gene type
expression_data <- data %>%
  filter(Type %in% c("Donor", "Native"), !is.na(deg_or_stab))
model <- lmer(logCount ~ Type * deg_or_stab + (1 | LGT.Cluster/LGT_native_accession), data = expression_data)
summary(model)

## make a figure 
ggplot(expression_data, aes(x = Type, y = logCount, fill = Type)) +
  geom_boxplot() +
  facet_wrap(~deg_or_stab) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Native" = "#76b0b0ff",
      "Donor"  = "#ab9fffff"
    )
  )
