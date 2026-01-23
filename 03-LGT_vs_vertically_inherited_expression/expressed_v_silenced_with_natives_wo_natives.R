library(dplyr)
library(tidyr)
library(textshape)
all_LGTs <- read.csv("02-LGT_expression/LGTs_only_silenced_expressed.csv") %>% 
  as_tibble() %>% 
  group_by(LGT.Cluster, Accession, Gene) %>% 
  summarise(Expression = first(Expression))
LGT_nat <- read.csv("01-Data/01-Lists/LGT & Natives.csv") %>% 
  as_tibble() %>% 
  filter(Type == "LGT") %>% 
  select(-Type, -Donor) %>% 
  mutate(with_native = "Yes")

df <- full_join(all_LGTs, LGT_nat) %>%
  mutate(with_native = ifelse(is.na(with_native), "No", with_native)) 

# do most LGTs have a native?
percen_with_native <- df %>% 
  group_by(Accession) %>% 
  summarise(
    total = n(),
    with_native = sum(with_native == "Yes")
  ) %>% 
  mutate(percentage_with_native = (with_native / total) * 100) %>% 
  ungroup()
percen_with_native %>% 
  summarise(
    mean_percen_w_native = mean(percentage_with_native),
    sd_percen_w_native = sd(percentage_with_native)
  )
  
# Create summary table
summary_table <- df %>%
  group_by(Accession, Expression) %>%
  summarise(
    total = n(),
    with_native_yes = sum(with_native == "Yes")
  ) %>%
  mutate(percentage_with_native = (with_native_yes / total) * 100) %>% 
  arrange(Expression, Accession)

# overall
summary_table %>% 
  group_by(Expression) %>% 
  summarise(mean_percen_w_native = mean(percentage_with_native),
            sd_percen_w_native = sd(percentage_with_native))

# stats test
# overall
# Create contingency table
contingency_table <- df %>%
  group_by(Expression, with_native) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = with_native, values_from = count, values_fill = 0) %>%
  column_to_rownames("Expression") # Convert first column to row names

# Perform Fisher's Exact Test
chi.test <- chisq.test(contingency_table)

# Print test result
print(chi.test)

