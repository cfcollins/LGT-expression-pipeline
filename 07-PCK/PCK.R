##### closer look at PCK 

## load libraries
library(ggplot2)
library(dplyr)

## load data
PCK <- read.csv("LGT_native_donor_counts_edit.csv") %>% 
  as_tibble() %>% 
  filter(LGT.Cluster == "LGT-019") %>% 
  mutate(Tissue_type = ifelse(Sample.organ == "Root", "Root", "Leaf/Sheath"))

# ZAM
ZAM <- PCK %>% 
  filter(LGT_native_accession == "ZAM") 
ZAM$Type <- factor(ZAM$Type, levels = c("LGT", "Native", "Donor"))
ZAM$Tissue_type <- factor(ZAM$Tissue_type, levels = c("Leaf/Sheath", "Root"))
ggplot(ZAM, aes(x = Type, y = logCount, shape = Tissue_type, colour = Gene)) +
  geom_jitter(width = .2, size = 2) +
  theme_minimal()

model <- lm(logCount ~ Type * Tissue_type, data = ZAM)
summary(model)

ZAM_LGT_roots <- ZAM %>% 
  filter(Type == "LGT" & Sample.organ == "Root")
ZAM_nat_don_roots <- ZAM %>% 
  filter(Type %in% c("Native", "Donor") & Sample.organ == "Root")
t.test(ZAM_LGT_roots$logCount, ZAM_nat_don_roots$logCount)

ZAM_nat_ls <- ZAM %>% 
  filter(Type == "Native" & Sample.organ != "Root")
ZAM_don_ls <- ZAM %>% 
  filter(Type == "Donor" & Sample.organ != "Root")
t.test(ZAM_nat_ls$logCount, ZAM_don_ls$logCount)

ZAM_nat_root <- ZAM %>% 
  filter(Type == "Native" & Sample.organ == "Root")
ZAM_don_root <- ZAM %>% 
  filter(Type == "Donor" & Sample.organ == "Root")
t.test(ZAM_nat_root$logCount, ZAM_don_root$logCount)

# AUS1
AUS1 <- PCK %>% 
  filter(LGT_native_accession == "AUS1") 
AUS1$Type <- factor(AUS1$Type, levels = c("LGT", "Native", "Donor"))
AUS1$Tissue_type <- factor(AUS1$Tissue_type, levels = c("Leaf/Sheath", "Root"))
ggplot(AUS1, aes(x = Type, y = logCount, fill = Tissue_type, colour = Gene)) +
  geom_boxplot() +
  theme_minimal()

model <- lm(logCount ~ Type * Tissue_type, data = AUS1)
summary(model)

AUS1_LGT_roots <- AUS1 %>% 
  filter(Type == "LGT" & Sample.organ == "Root")
AUS1_nat_don_roots <- AUS1 %>% 
  filter(Type %in% c("Native", "Donor") & Sample.organ == "Root")
t.test(AUS1_LGT_roots$logCount, AUS1_nat_don_roots$logCount)

AUS1_nat_ls <- AUS1 %>% 
  filter(Type == "Native" & Sample.organ != "Root")
AUS1_don_ls <- AUS1 %>% 
  filter(Type == "Donor" & Sample.organ != "Root")
t.test(AUS1_nat_ls$logCount, AUS1_don_ls$logCount)

AUS1_nat_root <- AUS1 %>% 
  filter(Type == "Native" & Sample.organ == "Root")
AUS1_don_root <- AUS1 %>% 
  filter(Type == "Donor" & Sample.organ == "Root")
t.test(AUS1_nat_root$logCount, AUS1_don_root$logCount)

# L04
L04 <- PCK %>% 
  filter(LGT_native_accession == "L04") 
L04$Type <- factor(L04$Type, levels = c("LGT", "Native", "Donor"))
L04$Tissue_type <- factor(L04$Tissue_type, levels = c("Leaf/Sheath", "Root"))
ggplot(L04, aes(x = Type, y = logCount, fill = Tissue_type, colour = Gene)) +
  geom_boxplot() +
  theme_minimal()

model <- lm(logCount ~ Type * Tissue_type, data = L04)
summary(model)

L04_LGT_nat <- L04 %>% 
  filter(Type %in% c("LGT", "Native"))
L04_don_ls <- L04 %>% 
  filter(Type == "Donor")
t.test(L04_LGT_nat$logCount, L04_don$logCount)

PCK$LGT_native_accession <- factor(PCK$LGT_native_accession, levels = c("AUS1", "ZAM", "L04"))
ggplot(PCK, aes(x = Type, y = logCount, shape = Tissue_type, colour = Accession)) +
  geom_boxplot(aes(fill = Tissue_type), alpha = .3, outlier.shape = NA) +
  geom_jitter(width = .2, size = 3, alpha = .8) +
  scale_fill_manual(values = c("white", "black")) +
  facet_wrap(~LGT_native_accession) +
  theme_bw()

ggsave("007-PCK/PCK.svg")
