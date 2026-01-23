### demonstrating PCK

## Load Libraries
library(ggplot2)
library(dplyr)

## Load Data
counts <- read.csv("COUNT_LISTS/LGT_natives_donors_counts.csv")
head(counts)

# filter down to PCK
PCK <- counts %>% 
  filter(LGT.Cluster == "LGT-019" & Type %in% c("LGT", "Native", "Donor"))

## Visualise
ggplot(PCK, aes(x = Type, y = logCount, fill = Accession)) +
  geom_boxplot()

ggplot(PCK, aes(x = Type, y = logCount, fill = Sample.organ)) +
  geom_boxplot()
