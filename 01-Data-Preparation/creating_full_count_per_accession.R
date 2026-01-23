### libraries
library(dplyr)
library(gplots)
library(ggplot2)
library(tidyr)

### load data

# count matrices
new_data <- read.csv("clean_salmon_matrices/new_data/align_to_own_ref/all_count_matrices/L04B_new_own_ref_TPM_count_matrix.csv")
old_data <- read.csv("clean_salmon_matrices/old_data/align_to_own_ref/all_count_matrices/L04B_old_own_ref_TPM_count_matrix.csv")
head(new_data)
head(old_data)

# sample info
sample_info <- read.csv("Expression project - Samples.csv")

### clean data

# merge dataframes 
combined_data <- merge(new_data, old_data, by = 'X', all = TRUE)
#combined_data <- new_data # this is for MRL
head(combined_data)

#result <- subset(combined_data, Gene == "ASEM_ZAM15-05-10_10714")

# need to remove the 'RA's (only for AUS1 and MRL)
combined_data$X <- base::sub("-RA$", "", combined_data$X)

# rename to 'X' column to 'Gene'
colnames(combined_data)[colnames(combined_data) == "X"] <- "Gene"
head(combined_data)

# change 'C4' to 'AUS1' in gene name (only for AUS1)
# this needs changing for different accessions

combined_data$Gene <- gsub("ASEM_C4_", "ASEM_AUS1_", combined_data$Gene)
head(combined_data)

# transpose the data 
transposed_data <- pivot_longer(combined_data, cols = -Gene, names_to = "Sample", values_to = "Count")
all_data <- merge(transposed_data, sample_info, by = 'Sample', all = FALSE)

### make spreadsheet

write.csv(all_data, "clean_salmon_matrices/full_accession_matrices/L04B_full_count_matrix.csv")
