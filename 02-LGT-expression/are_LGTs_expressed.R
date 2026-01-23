#### ARE LGTs EXPRESSED? ####
#### CAT COLLINS ############
#### REWRITTEN APR 2024 #####

## Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(svglite)
library(ggridges)

## Load data
LGTs <- read.csv("001-Data/002-Counts/LGTs_only_counts.csv") %>% 
  as_tibble() %>% 
  select(!c(Representative, Identification, Fragment))

## Labeling for silenced / expressed ####

# we're going to use the cut off as less than 0.5 TPM for silenced
# if all the samples for a gene is less than 0.5, then that gene is 
# categorised as 'silenced'
# if any of the samples are 0.5 or over, then that gene is categorised as 
# 'expressed'

# expressed or silenced genes
silenced_expressed <- LGTs %>% 
  group_by(Gene) %>% 
  mutate(Expression = ifelse(any(Count >= 0.5), "Expressed", "Silenced")) %>% 
  ungroup()
write.csv(silenced_expressed, "02-LGT_expression/LGTs_only_silenced_expressed.csv", row.names = FALSE)

# Step 1: Calculate mean and SD of logCount for each Gene
gene_stats <- silenced_expressed %>%
  group_by(Gene) %>%
  summarize(
    GeneMeanLogCount = mean(logCount, na.rm = TRUE),
    GeneSDLogCount = sd(logCount, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Join the per-gene stats back to the original dataframe
df_with_gene_stats <- silenced_expressed %>%
  left_join(gene_stats, by = "Gene")

# Step 3: Calculate summary statistics for each Accession
unique_genes_per_accession <- df_with_gene_stats %>%
  group_by(Accession) %>%
  summarize(
    UniqueGeneCount = n_distinct(Gene),
    ExpressedGeneCount = n_distinct(Gene[Expression == "Expressed"]),
    ExpressedPercentage = (ExpressedGeneCount / UniqueGeneCount) * 100,
    MeanLogCount = mean(GeneMeanLogCount, na.rm = TRUE),
    MinLogCount = min(logCount),
    maxLogCount = max(logCount),
    SDLogCount = sd(GeneMeanLogCount, na.rm = TRUE)  # SD across gene means
  )

# mean and standard deviation for expression percentages
data <- c(80.8, 75, 76.5, 88.9, 66.1)
sd(data)
mean(data)

write.csv(unique_genes_per_accession, "02-LGT_expression/summary_stats_per_accession.csv", row.names = FALSE)

## visualisation on a point plot #####

# reorder Accession
LGTs$Accession = factor(LGTs$Accession, levels = c("AUS1", "ZAM", "L04", "KWT", "MRL"))

# for each sample within each accession, this counts the no. of genes with
# a count over 0.5, then divides it by the total no. of genes to get a %
percentages_by_sample <- LGTs %>% 
  group_by(Accession, Sample.organ, Sample) %>% 
  summarise(percentage_expressed = sum(Count >= 0.5) / n() * 100) %>% 
  ungroup()

percentages_by_sample %>% 
  group_by(Accession) %>% 
  summarise(mean = mean(percentage_expressed))
# this checks for each gene in each accession whether any of the samples have
# a count that is over 0.5 TPM. If any of them have this, it is marked 
# 'Expressed'
# if all samples for this gene have a count less than 0.5, it is marked 
# 'Silenced'
# this then adds together the genes that are expressed for each accession
# and uses the total to create a percentage
overall <- LGTs %>% 
  group_by(Gene) %>% 
  mutate(Expression = ifelse(any(Count >= 0.5), "Expressed", "Silenced")) %>% 
  group_by(Accession) %>% 
  summarise(percentage = round(sum(Expression == "Expressed") / n() * 100, 1)) %>% 
  mutate(text = paste(percentage, "%", sep = "")) %>% 
  ungroup() %>%
  mutate(Accession = recode(Accession, 
                            "MRL" = "AANG_UGA4", 
                            "ZAM" = "ZAM1505-10", 
                            "KWT" = "RSA5-3", 
                            "L04" = "TAN1-04B",
                            "AUS1" = "AUS1"))

# change the accession names
percentages_by_sample <- percentages_by_sample %>%
  mutate(Accession = recode(Accession, 
                            "MRL" = "AANG_UGA4", 
                            "ZAM" = "ZAM1505-10", 
                            "KWT" = "RSA5-3", 
                            "L04" = "TAN1-04B",
                            "AUS1" = "AUS1"))

percentages_by_sample %>% 
  distinct(Accession)
# this creates a boxplot out of the percentages per sample per accession
# it then overlays points for sample.organ using the percentage by organ table
# it then overlays points for the percentage of expressed genes across each 
# accession with crosses on them
# it then adds labels to each of these points and adds a red line at 50
ggplot(percentages_by_sample, aes(x = percentage_expressed, y = factor(Accession), fill = Accession)) +
  geom_boxplot() +
  geom_point(data = overall, aes(x = percentage, y = Accession, colour = Accession), size = 12) +
  geom_point(data = overall, aes(x = percentage, y = Accession), colour = "black", size = 12, pch = 4, show.legend = FALSE) +
  geom_text(data = overall, aes(x = percentage, y = Accession, label = text), vjust = -1, size = 8) +  # Increase text size
  geom_vline(xintercept = 50, colour = "red", linetype = "dashed", size = 1) +
  guides(fill = "none", colour = "none") +  # Remove legends
  xlim(0, 100) +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 5)) +
  scale_colour_manual(values = wes_palette("Darjeeling1", n = 5)) +
  theme_ridges(font_size = 20) +  # Increase base font size
  theme(axis.ticks.y = element_blank()) + # Remove y-axis ticks
  labs(x = "Percentage Expressed")


ggplot(percentages_by_sample, aes(x = percentage_expressed, y = factor(Accession), fill = Accession)) +
  geom_boxplot() +
  geom_point(data = overall, aes(x = percentage, y = Accession, colour = Accession), size = 12) +
  geom_point(data = overall, aes(x = percentage, y = Accession), colour = "black", size = 13, pch = 4, show.legend = FALSE) +
  geom_vline(xintercept = 50, colour = "red", linetype = "dashed", size = 1) +
  guides(fill = "none", colour = "none") +  # Remove legends
  xlim(0, 100) +
  theme_minimal() +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 5)) +
  scale_colour_manual(values = wes_palette("Darjeeling1", n = 5)) +  # Increase base font size
  theme(axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank(),   # Remove y-axis text
        axis.ticks.y = element_blank()) + # Remove y-axis ticks
  labs(x = "Percentage Expressed")



ggsave("002-LGT-expression/LGT_expression_boxplots.svg", width = 10, height = 6.99, dpi = 1200)


ggsave("2. Are LGTs Expressed/final_plot_w_histograms.png", combined_plot1, width = 12, height = 7)
combined_plot2 <- (plot)
ggsave("2. Are LGTs Expressed/final_plot.png", combined_plot2, dpi = 800, width = 10, height = 4)
ggsave("02-LGT_expression//final_plot.svg", combined_plot2, dpi = 800, width = 10, height = 4)

svglite("2. Are LGTs Expressed/final_plot.svg", width = 8, height = 6)
print(combined_plot2)  # Print the plot to save it
dev.off() 

svglite("2. Are LGTs Expressed/final_plot_part1.svg", width = 8, height = 6)
print(tree_plot)  # Print the plot to save it
dev.off() 
svglite("2. Are LGTs Expressed/final_plot_part2.svg", width = 8, height = 6)
print(plot)  # Print the plot to save it
dev.off() 
