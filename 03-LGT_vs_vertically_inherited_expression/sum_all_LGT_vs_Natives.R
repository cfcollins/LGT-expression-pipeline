library(dplyr)
library(ggplot2)
library(wesanderson)
library(broom)
library(svglite)

# accessing data 
counts <- read.csv("001-Data/002-Counts/LGT_natives_counts.csv") %>% 
  filter(Type == "LGT" | Type == "Native") %>% 
  as_tibble() %>% 
  group_by(LGT.Cluster, Type, Sample.organ, Sample, Accession) %>% 
  summarise(logCount = sum(logCount))

counts$Accession = factor(counts$Accession, levels = c("MRL", "KWT", "L04", "AUS1", "ZAM"))

# creating differences
counts_edit <- counts %>% 
  left_join(counts, by = c("LGT.Cluster", "Sample", "Accession")) %>% 
  filter(Type.x == "LGT" & Type.y == "Native") %>% 
  mutate(logCount_diff = logCount.x - logCount.y) %>% 
  select(-logCount.x, -logCount.y)

# visualising 
ggplot(counts_edit, aes(x = logCount_diff, y =  reorder(LGT.Cluster, logCount_diff), fill = Accession)) +
  geom_boxplot(linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5)) +
  labs(x = "logCount Native subtracted from logCount LGT",
       y = "LGT Cluster") +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 5)) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# stats tests
counts <- group_by(counts, LGT.Cluster, Accession)
results <- do(counts, tidy(wilcox.test(.$logCount ~ .$Type,
                                       alternative = "two.sided")))

p_values <- results$p.value
corrected_p_values <- p.adjust(p_values, method = "BH")
results$corrected_p_value <- corrected_p_values
results <- results %>% 
  mutate(Significance = ifelse(is.na(corrected_p_value), "Non-Significant", 
                               ifelse(corrected_p_value < 0.05, "Significant", "Non-Significant")))
write.csv(results, "003-LGT_vs_native_expression/02-sum_all/t-test_results_diff_sum_all.csv", row.names = FALSE)
sig <- results %>% 
  filter(Significance == "Significant")
nonsig <- results %>% 
  filter(Significance == "Non-Significant")

sig_counts <- semi_join(counts_edit, sig, by = c("LGT.Cluster", "Accession")) %>% 
  mutate(Significance = "Significant")
nonsig_counts <- semi_join(counts_edit, nonsig, by = c("LGT.Cluster", "Accession")) %>% 
  mutate(Significance = "Non-Significant")

counts_w_stats <- rbind(sig_counts, nonsig_counts)

edit <- counts_w_stats %>% 
  group_by(LGT.Cluster, Accession) %>% 
  summarise(Significance = first(Significance),
            logCount_diff = mean(logCount_diff)) %>% 
  mutate(group = case_when(
    Significance == "Significant" & logCount_diff < 0 ~ "Significantly less",
    Significance == "Significant" & logCount_diff > 0 ~ "Significantly more",
    Significance == "Non-Significant" ~ "No difference"
  ))

summary <- edit %>% 
  group_by(Accession) %>% 
  summarise(
    total = n(),
    sig_less_n = sum(group == "Significantly less"),
    percen_sig_less = (sig_less_n / total) * 100 
  )
write.csv(summary, "003-LGT_vs_native_expression/02-sum_all/summary_genes_LGT_vs_Natives.csv", row.names = FALSE)

edit2 <- edit %>% 
  group_by(LGT.Cluster) %>% 
  summarise(Consistency = case_when(
    all(group == "Significantly less") ~ "Significantly less",
    all(group == "Significantly more") ~ "Significantly more",
    all(group == "No difference") ~ "No difference",
    TRUE ~ "Inconsistent"
  )) %>% 
  ungroup() %>% 
  count(Consistency)

# visualise again

counts_w_stats <- counts_w_stats %>%
  mutate(Significance = factor(Significance, levels = c("Significant", "Non-Significant")))

# let's try combining the LGT Cluster and Accession
counts_combined_columns <- counts_w_stats %>%
  mutate(combined_column = paste(LGT.Cluster, Accession, sep = "_"))
colours <- c("Significant" = "#72bf6a",
             "Non-Significant" = "#ff2c2c")

counts_combined_columns %>% 
  group_by(LGT.Cluster, Accession) %>% 
  filter(all(logCount_diff == 0)) %>% 
  print(n = 33)

ggplot(counts_combined_columns %>% filter(!LGT.Cluster %in% c("LGT-004", "LGT-039", "LGT-174")), aes(x = logCount_diff, 
                                    y = reorder(combined_column, logCount_diff), 
                                    fill = Significance)) +
  geom_boxplot(linewidth = 0.2, outlier.size = 0.5) +  # Reduce outlier point size
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
  scale_fill_manual(values = colours) + 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_line(size = 0.2),  # Reduce major grid line width
    panel.grid.minor = element_line(size = 0.1),  # Reduce minor grid line width
    legend.position = "none"  # Remove legend
  )


ggsave("003-LGT_vs_native_expression/02-sum_all/LGT_vs_natives_sum_all.svg", width = 6, height = 5)

write.csv(counts_combined_columns, "003-LGT_vs_native_expression/02-sum_all/sum_all_sig_marked.csv", row.names = FALSE)

## highlighting PCK and PEPC
editing <- counts_combined_columns %>% 
  mutate(note = case_when(
    LGT.Cluster == "LGT-019" ~ "PCK",
    LGT.Cluster == "LGT-079" ~ "PEPC",
    TRUE ~ "none"
  ))

ggplot(editing, aes(x = logCount_diff,
                    y =  reorder(combined_column, logCount_diff),
                    fill = note)) +
  geom_boxplot(linewidth = 0.3) +
  theme(axis.text.y = element_text(size = 4))

## STATS

### linear mixed effects model
library(lme4)
library(lmerTest)

counts_combined_columns

model <- lmer(logCount_diff ~ 1 + (1 | LGT.Cluster / Accession), data = counts_combined_columns)
summary(model)

t_value <- fixef(model)["(Intercept)"] / summary(model)$coefficients[, "Std. Error"]
p_value <- pt(t_value, df.residual(model), lower.tail = TRUE)
print(p_value)

## plot for figure
summary <- counts_combined_columns %>% 
  group_by(combined_column) %>% 
  summarise(median = median(logCount_diff),
            lower = quantile(logCount_diff, 0.25, na.rm = TRUE),
            upper = quantile(logCount_diff, 0.75, na.rm = TRUE)) %>% 
  filter(!(median == 0 & lower == 0 & upper == 0))


ggplot(summary, aes(x = median, y = reorder(combined_column, median))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 1) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, color = "darkgrey", linewidth = 0.5) +
  geom_point(size = .8, color = "black") +
  theme(
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.y = element_blank()   
  )

ggsave("003-LGT_vs_native_expression/02-sum_all/figure1b.svg", width = 6, height = 6)

# "LGT-004", "LGT-039", "LGT-174"
counts_combined_columns %>% 
  filter(LGT.Cluster == "LGT-039") %>% 
  print(n = 35)

counts_combined_columns %>% 
  group_by(combined_column) %>% 
  summarise(median = median(logCount_diff),
            lower = quantile(logCount_diff, 0.25, na.rm = TRUE),
            upper = quantile(logCount_diff, 0.75, na.rm = TRUE)) %>% 
  filter(median == 0 & lower == 0 & upper == 0)

counts_combined_columns %>% 
  group_by(combined_column) %>% 
  summarise(median = median(logCount_diff),
            lower = quantile(logCount_diff, 0.25, na.rm = TRUE),
            upper = quantile(logCount_diff, 0.75, na.rm = TRUE)) %>% 
  filter(median == 0)
