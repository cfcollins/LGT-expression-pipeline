# removing LGTs where both homologs are silenced

## load libraries
library(dplyr)

## load data
all_counts <- read.csv("LGT_native_donor_counts_edit.csv") %>% 
  as_tibble() 
homologs_silenced <- all_counts %>% 
  filter(Type %in% c("Native", "Donor")) %>% 
  group_by(LGT.Cluster, LGT_native_accession) %>% 
  mutate(Expression = ifelse(all(Count < 0.5), "Silenced", "Expressed")) %>% 
  filter(Expression == "Silenced") %>% 
  group_by(LGT.Cluster, Type) %>% 
  distinct(Gene) %>% 
  arrange(LGT.Cluster) %>% 
  print(n = 21)

# LGT-039
# LGT-078
# LGT-142
# LGT-177

counts_filtered <- all_counts %>% 
  filter(!LGT.Cluster %in% c("LGT-039", "LGT-078", "LGT-142", "LGT-177"))

write.csv(counts_filtered, "LGT_native_donor_counts_filtered.csv", row.names = FALSE)

# how many genes have active and expressed native and donor homologs
all_LGTs <- read.csv("01-Data/01-Lists/LGTs only.csv") %>% 
  select(-Representative, -Identification, -Fragment)
combined <- counts_filtered %>% 
  filter(Type == "LGT") %>% 
  group_by(LGT.Cluster, Accession, Gene) %>% 
  summarise(logCount = mean(logCount)) %>% 
  mutate(A_E_nat_don = "Yes") %>% 
  full_join(all_LGTs, by = c("LGT.Cluster", "Accession", "Gene")) %>% 
  mutate(A_E_nat_don = ifelse(is.na(A_E_nat_don), "No", A_E_nat_don)) 

# by accession - clusters
combined %>% 
  group_by(LGT.Cluster, Accession) %>% 
  summarise(A_E_nat_don = first(A_E_nat_don)) %>% 
  group_by(Accession) %>% 
  summarise(
    n = n(),
    A_E_nat_don_yes = sum(A_E_nat_don == "Yes"),
    percen_yes = (A_E_nat_don_yes / n)*100
  )

# by accession - including duplicates
combined %>% 
  group_by(Accession) %>% 
  summarise(
    n = n(),
    A_E_nat_don_yes = sum(A_E_nat_don == "Yes"),
    percen_yes = (A_E_nat_don_yes / n)*100
  )

# across all LGTs 
combined %>% 
  group_by(LGT.Cluster) %>% 
  summarise(A_E_nat_don = ifelse(any(A_E_nat_don == "Yes"), "Yes", "No")) %>% 
  count(A_E_nat_don)
