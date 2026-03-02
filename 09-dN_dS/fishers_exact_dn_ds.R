## fisher's exact test for dN/dS between LGT and vertically inherited genes

# libraries
library(dplyr)
library(ggplot2)

# data
LGT_dn_ds <- read.csv("LGT_natives_dn_ds.csv") %>% 
  as_tibble() %>% 
  filter(Type == "LGT")
categories <- read.csv("expanded_deg_or_stab_SUM_ALL.csv")

native_dn_ds <- read.csv("LGT_natives_dn_ds.csv") %>% 
  as_tibble() %>% 
  filter(Type == "Native")

LGT_dn_ds %>% 
  count(Classification)
native_dn_ds %>% 
  count(Classification)

# Make the 2×2 table
tbl <- matrix(
  c(19, 97,   # LGT: relaxed, not relaxed
    2, 114),  # Native: relaxed, not relaxed
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Category = c("LGT", "Native"),
    Constraint = c("Relaxed", "Not_relaxed")
  )
)

# Fisher's exact test
fisher.test(tbl)

### degenerating vs put stab
LGT_dn_ds %>% 
  full_join(categories) %>%
  filter(!is.na(Classification),
         !is.na(deg_or_stab)) %>% 
  group_by(deg_or_stab) %>% 
  count(Classification)

# Make the 2×2 table
tbl <- matrix(
  c(17, 66,   # deg: relaxed, not relaxed
    8, 44),  # put stab: relaxed, not relaxed
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Category = c("deg", "put_stab"),
    Constraint = c("Relaxed", "Not_relaxed")
  )
)

# Fisher's exact test
fisher.test(tbl)
