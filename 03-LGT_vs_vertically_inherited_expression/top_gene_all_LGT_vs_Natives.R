library(dplyr)
library(ggplot2)
library(wesanderson)
library(broom)
library(svglite)

# accessing data 
counts <- read.csv("001-Data/002-Counts/LGT_natives_counts.csv") %>% 
  as_tibble() %>% 
  filter(Type == "LGT" | Type == "Native") %>% 
  group_by(LGT.Cluster, Type, Sample, Accession) %>% 
  summarise(logCount = max(logCount))


setwd("003-LGT_vs_native_expression/03-top_gene_all/")
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
write.csv(results, "t-test_results_diff_top_gene_all.csv", row.names = FALSE)
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
write.csv(summary, "summary_genes_LGT_vs_Natives.csv", row.names = FALSE)

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

ggplot(counts_combined_columns, aes(x =logCount_diff, 
                           y =  reorder(combined_column, logCount_diff), 
                           fill = Significance)) +
  geom_boxplot(linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 2) + 
  scale_fill_manual(values = colours) + 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank() 
  )
ggsave("LGT_vs_natives_top_gene_all.png", dpi = 1200)

write.csv(counts_combined_columns, "top_gene_all_sig_marked.csv", row.names = FALSE)

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
