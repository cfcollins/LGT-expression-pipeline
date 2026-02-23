## script to look at whether truncated vertically inherited genes have
## full length LGTs

# libraries
library(dplyr)
library(tidyr)

# data
trunc <- read.csv("v2_pseudogene_LGT_vs_native_adjusted_for_native_reuse.csv") %>% 
  as_tibble()
data <- trunc %>% 
  group_by(LGT.Cluster, LGT_native_accession, Type) %>% 
  summarise(Truncation = ifelse(all(Pseudogene == "Yes"), "all_truncated", "fine")) %>% 
  pivot_wider(id_cols = c("LGT.Cluster", "LGT_native_accession"), 
              names_from = Type, 
              values_from = Truncation) %>% 
  summarise(Truncation = case_when(
    LGT == "fine" & Native == "fine" ~ "all_fine",
    LGT == "all_truncated" & Native == "fine" ~ "LGT truncated, native fine",
    LGT == "fine" & Native == "all_truncated" ~ "LGT fine, native truncated",
    LGT == "all_truncated" & Native == "all_truncated" ~ "all truncated"
  )) %>% 
  ungroup()

data %>% 
  count(Truncation)

            
# I want to compare this to the overall LGT truncated rate 
binom.test(x = 10, n = 17, p = 0.667)

