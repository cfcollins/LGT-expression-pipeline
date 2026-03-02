### dN comparison
data <- read.table("011-divergence_native_vs_donor/02-genetic divergence/results.txt") %>% 
  as_tibble() %>% 
  rename(LGT.Cluster = "V1",
         Accession = "V2",
         S = "V5",
         N = "V6",
         t = "V7",
         kappa = "V8",
         omega = "V9",
         dN = "V10",
         dN_SE = "V11",
         dS = "V12",
         dS_SE = "V13")
ggplot(data, aes(x = log(log(dS+1)))) +
  geom_histogram()
ggplot(data, aes(x = log(log(dN+1)))) +
  geom_histogram()
ggplot(data, aes(x = log(log(omega+1)))) +
  geom_histogram()


write.csv(data, "genetic_divergence_native_donor.csv", row.names = FALSE)

categories <- read.csv("expanded_deg_or_stab_SUM_ALL.csv") %>% 
  rename(Accession = "LGT_native_accession")

categories

combined <- full_join(data,categories, by = c("LGT.Cluster", "Accession")) %>% 
  filter(!is.na(deg_or_stab))

## average dS for degenerating and put stab
combined %>% 
  group_by(deg_or_stab) %>% 
  summarise(mean_dS = mean(dS),
            SD_dS = sd(dS))

combined %>% 
  select(LGT.Cluster, Accession, deg_or_stab, dS) %>% 
  print(n = 148)

ggplot(combined, aes(x = deg_or_stab, y = log(log(dS + 1)), fill = deg_or_stab)) +
  geom_boxplot() +
  facet_wrap(~Accession) +
  theme_minimal()

log(dN + 1)
log(log(dN + 1))


model <- lm(dS ~ deg_or_stab, data = combined)
summary(model)
model <- lmer(dS ~ deg_or_stab + (1 | Accession), data = combined)
summary(model)

wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = combined)
wilcox.test(dN ~ deg_or_stab, data = combined)
wilcox.test(omega ~ deg_or_stab, data = combined)

dS_model <- lmer(log(log(dS+1)) ~ deg_or_stab + (1 | Accession), data = combined)
summary(dS_model)
dN_model <- lmer(log(log(dN+1)) ~ deg_or_stab + (1 | Accession), data = combined)
summary(dN_model)
omega_model <- lmer(log(log(omega+1))  ~ deg_or_stab + (1 | Accession), data = combined)
summary(omega_model)

log(dN+1)
log(log(dS+1))
ggplot(combined, aes(x = deg_or_stab, y = log(log(dS+1)), fill = deg_or_stab)) +
  geom_boxplot() +
  theme_minimal()
wilcox.test(dS ~ deg_or_stab, data = combined)

ggplot(combined %>% filter(dS < 1), aes(x = deg_or_stab, y = dS, fill = deg_or_stab)) +
  geom_boxplot()



###
exp_cat <- read.csv("LGTs_less_than_native_and_donor_or_not-SUM-ALL.csv") %>% 
  rename(Accession = "LGT_native_accession")
combined <- full_join(data, exp_cat, by = c("LGT.Cluster", "Accession")) %>% 
  filter(!is.na(LGT_lower_than_homologs))


ggplot(combined, aes(x = LGT_lower_than_homologs, y = log(log(dS+1)), fill = LGT_lower_than_homologs)) +
  geom_boxplot() +
  facet_wrap(~Accession)


model <- lm(log(log(dS+1)) ~ LGT_lower_than_homologs, data = combined)
summary(model)

library(lme4)
library(lmerTest)
model <- lm(log(log(dS+1)) ~ LGT_lower_than_homologs * Accession, data = combined)
summary(model)

########

AUS1 <- combined %>% 
  filter(Accession == "AUS1")
KWT <- combined %>% 
  filter(Accession == "KWT")
L04 <- combined %>% 
  filter(Accession == "L04")
ZAM <- combined %>% 
  filter(Accession == "ZAM")
MRL <- combined %>% 
  filter(Accession == "MRL")

wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = AUS1)
wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = KWT)
wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = L04)
wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = ZAM)
wilcox.test(log(log(dS+1)) ~ deg_or_stab, data = MRL)


ggplot(combined, aes(x = deg_or_stab, y = log(log(dS + 1)), fill = deg_or_stab)) +
  geom_boxplot() +
  facet_wrap(~Accession) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Putatively stable" = "#1a9850ff",
      "Degenerating"  = "#f46d43ff"
    )
  )
ggsave("genetic_divergence.svg")
