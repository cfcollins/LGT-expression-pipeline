#########################################################################################
#       Script Name: pca.R
#       Description: PCA for rna-seq data aligned to A. semialata genome
#       Author:      Ben Alston
#       Date:        March 2023
#########################################################################################

library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggtext)
library(rtracklayer)
library(ggpubr)
library(dplyr)

#load raw data ----

new_count_matrix <- read.csv('1. Data Preparation/1. RNA-Seq/4. Sensible dataframe generation in R/
                             clean_salmon_matrices/new_data/align_to_AUS1/all_new_AUS1_TPM_count_matrix.csv')
head(new_count_matrix)
old_count_matrix <- read.csv('clean_salmon_matrices/old_data/align_to_AUS1/all_old_AUS1_TPM_count_matrix.csv')
head(old_count_matrix)
# extract colnames to change count matrix colnames
# new
all_sampleinfo1 <- read.csv('Bens stuff/LGT_count_lists/all_sampleinfo.csv', row.names = NULL)
head(all_sampleinfo)
all_sampleinfo <- read.csv("Expression Project - Sample Overview - Samples.csv", row.names = NULL) %>% 
  as_tibble() %>% 
  select(Data, Dataset, Species, Accession, Sample.organ)

all_matrix <- cbind(new_count_matrix, old_count_matrix) %>% 
  column_to_rownames(var = 'X') %>% 
  select(-X) %>% 
  as.matrix()


dds_all <- DESeqDataSetFromMatrix(round(all_matrix), colData = all_sampleinfo, design = ~Dataset)
dds_all <- estimateSizeFactors(dds_all)
vst_all <- vst(dds_all, blind = TRUE)
raw_PCA <- plotPCA(vst_all, intgroup = 'Dataset')

plotPCA(vst_all, intgroup = 'Dataset')

initial_PCA_data <- plotPCA(vst_all, intgroup = 'Dataset', returnData = TRUE)
initial_PCA_data  <- initial_PCA_data  %>% select(-group, -Dataset)
initial_PCA_data <- right_join(initial_PCA_data, all_sampleinfo,  by = c('name' = "Data"))

new_only <- initial_PCA_data %>% 
  filter(Dataset == "New")
plot <- ggplot(initial_PCA_data, aes(x = PC1, y = PC2, col = Sample.organ, shape = Species)) +
  theme_bw() +
  geom_point(size = 3) +
  geom_point(data = new_only, aes(x = PC1, y = PC2), shape = 1, size = 7, colour = "grey") +
  labs(x = 'PC1: 44% variance', y = 'PC2: 19% variance',
       color = "Tissue Type")

ggsave("PCA_all.png", plot)

ggplot(initial_PCA_data, aes(x = PC1, y = PC2, col = Sample.organ, shape = Species)) +
  theme_bw() +
  geom_point(size = 3) + 
  labs(x = 'PC1: 82% variance', y = 'PC2: 7% variance',
       col = 'Accession', shape = 'Species') +  
  scale_colour_manual(values = c('#34d2eb', '#1b1bab', '#8108c2','#09fa05', '#b0e010', '#e01017', '#346beb')) +
  #scale_colour_manual(values = pokepal(pokemon ='pikachu', spread = 7)) +
  scale_shape_discrete(labels = c('*A. angusta*','*A. semialata*','*S. italica*', '*T. triandra*')) +
  theme(legend.text = element_markdown())

jpeg('plots/allpca2.jpg', units = 'cm', width = 30, height = 15, res = 600)
ggplot(initial_PCA_data, aes(x = PC1, y = PC2, col = accession, shape = species)) +
  theme_bw() +
  geom_point(size = 3) + 
  labs(x = 'PC1: 82% variance', y = 'PC2: 7% variance',
       col = 'Accession', shape = 'Species') +  
  scale_colour_manual(values = c('#34d2eb', '#1b1bab', '#8108c2','#09fa05', '#b0e010', '#e01017', '#346beb')) +
  #scale_colour_manual(values = pokepal(pokemon ='pikachu', spread = 7)) +
  scale_shape_discrete(labels = c('*A. angusta*','*A. semialata*','*S. italica*', '*T. triandra*')) +
  theme(legend.text = element_markdown())
dev.off()


sym_diff2 <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))
sym_diff2(colnames(all_matrix),
         all_sampleinfo$sample_id)



## just Alloteropsis semialata



just_sem <- all_sampleinfo %>% 
  filter(Species == "Alloteropsis semialata")

a_new <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/AUS1_new_AUS1_TPM_count_matrix.csv")
k_new <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/KWT_new_AUS1_TPM_count_matrix.csv")
l_new <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/L04B_new_AUS1_TPM_count_matrix.csv")
z_new <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/ZAM_new_AUS1_TPM_count_matrix.csv")
a_old <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/AUS1_old_AUS1_TPM_count_matrix.csv")
k_old <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/KWT_old_AUS1_TPM_count_matrix.csv")
l_old <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/L04B_old_AUS1_TPM_count_matrix.csv")
z_old <- read.csv("clean_salmon_matrices/just_sem_align_to_AUS1/ZAM_old_AUS1_TPM_count_matrix.csv")

matrix_sem <- cbind(a_new, k_new, l_new, z_new, a_old, k_old, l_old, z_old) %>% 
  column_to_rownames(var = 'X') %>% 
  select(-X) %>% 
  as.matrix()

dds_all <- DESeqDataSetFromMatrix(round(matrix_sem), colData = just_sem, design = ~Dataset)
dds_all <- estimateSizeFactors(dds_all)
vst_all <- vst(dds_all, blind = TRUE)
raw_PCA <- plotPCA(vst_all, intgroup = 'Dataset')

plotPCA(vst_all, intgroup = 'Dataset')

initial_PCA_data <- plotPCA(vst_all, intgroup = 'Dataset', returnData = TRUE)
initial_PCA_data  <- initial_PCA_data  %>% select(-group, -Dataset)
initial_PCA_data <- right_join(initial_PCA_data, just_sem,  by = c('name' = "Data"))

new_only <- initial_PCA_data %>% 
  filter(Dataset == "New")
plot <- ggplot(initial_PCA_data, aes(x = PC1, y = PC2, col = Accession, shape = Sample.organ)) +
  theme_bw() +
  geom_point(size = 3) +
  geom_point(data = new_only, aes(x = PC1, y = PC2), shape = 1, size = 7, colour = "grey") +
  labs(x = 'PC1: 55% variance', y = 'PC2: 13% variance',
       shape = "Tissue Type")

ggsave("PCA_sem_only.png", plot)
