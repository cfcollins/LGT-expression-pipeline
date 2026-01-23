# comparing whether a different proportion of genes are expressed for 
# LGTs with vertically inherited orthologs and
# LGTs without vertically inherited orthologs

## load libraries
library(dplyr)
library(ggplot2)

## load data
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

# stats
# per accession
per_accession <- df %>% 
  group_by(Accession, with_native) %>% 
  summarise(n = n(),
    Expressed = sum(Expression == "Expressed"),
    percen_expressed = Expressed / n)





chi_results <- per_accession %>%
  group_by(Accession) %>%
  summarise(
    chisq_p_value = {
      # Create 2x2 table: rows = with_native, columns = Expressed/Not expressed
      tab <- matrix(
        c(
          Expressed[with_native == "Yes"],
          n[with_native == "Yes"] - Expressed[with_native == "Yes"],
          Expressed[with_native == "No"],
          n[with_native == "No"] - Expressed[with_native == "No"]
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(with_native = c("Yes", "No"), Expression = c("Expressed", "NotExpressed"))
      )
      if (any(tab < 5)) {
        # Use Fisher's exact test if expected counts are small
        fisher.test(tab)$p.value
      } else {
        chisq.test(tab)$p.value
      }
    }
  )

# overall
data <- df %>% 
  group_by(with_native) %>% 
  summarise(n = n(),
            Expressed = sum(Expression == "Expressed"),
            percen_expressed = Expressed / n)
data %>% 
  summarise(
    chisq_p_value = {
      # Create 2x2 table: rows = with_native, columns = Expressed/Not expressed
      tab <- matrix(
        c(
          Expressed[with_native == "Yes"],
          n[with_native == "Yes"] - Expressed[with_native == "Yes"],
          Expressed[with_native == "No"],
          n[with_native == "No"] - Expressed[with_native == "No"]
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(with_native = c("Yes", "No"), Expression = c("Expressed", "NotExpressed"))
      )
      if (any(tab < 5)) {
        # Use Fisher's exact test if expected counts are small
        fisher.test(tab)$p.value
      } else {
        chisq.test(tab)$p.value
      }
    }
  )
